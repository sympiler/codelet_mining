#include "SPMV_CVR.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#define MICRO_IN_SEC 1000000.00
#include <omp.h>

#include <unistd.h>
#include <sys/mman.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <list>
#include <vector>
#include <set>
#define WRITE

#include <immintrin.h>

using namespace std;

#ifndef TYPE_PRECISION
#define TYPE_PRECISION

//#define TYPE_SINGLE
#define TYPE_DOUBLE

#define FLOAT_TYPE double

#ifdef TYPE_SINGLE
#define SIMD_LEN 16
#else
#define SIMD_LEN 8
#endif

#endif


#define FIELD_LENGTH 128
#define floatType double
//#define DONT_USE_MM_MALLOC
#define USE_MM_MALLOC
#define ALIGN 64
#ifdef USE_MM_MALLOC
#define ALLOC(t,s) (t *)_mm_malloc((s)*sizeof(t), ALIGN)
#define FREE(p) _mm_free(p)
#else
#define ALLOC(t,s) new t[s]
#define FREE(p) delete[] p
#endif
#ifdef __AVX512__

__m512i pbadc = _mm512_set_epi32(4,5,6,7,0,1,2,3,12,13,14,15,8,9,10,11);
#define _MM_PERM_BADC pbadc
#define _mm512_permute4f128_epi32(idx, p) _mm512_permutexvar_epi32(p, idx)

#define _MM_SCALE_4 4
#define _MM_SCALE_8 8

#define _mm512_i32logather_pd(vi, ba, s) _mm512_i32gather_pd(vi, ba, s)

namespace SPMV_CVR {
    struct Coordinate {
        int x;
        int y;
        float val;
    };


    class i_Range {

    public:
        int i_start;
        int i_end;

        i_Range(int a, int b) {
            i_start = a;
            i_end = b;
        }
        i_Range() {
            i_start = 0;
            i_end = 0;
        }
    };

    struct v_Range {
        int v_start;
        int v_end;
    };

    typedef __attribute__((aligned(64))) union zmmd {
        __m512d reg;
        __m512i regi32;
        double elems[8];
        int elemsi32[16];
    } zmmd_t;

    double microtime() {
        int tv_sec, tv_usec;
        double time;
        struct timeval tv;
        struct timezone tz;
        gettimeofday(&tv, &tz);

        return tv.tv_sec + tv.tv_usec / MICRO_IN_SEC;
    }

    class Edge {

    public:
        int el;
        int er;
        FLOAT_TYPE value;

        Edge(int a, int b, FLOAT_TYPE v) {
            el = a;
            er = b;
            value = v;
        }
        Edge(int a, int b) {
            el = a;
            er = b;
        }
    };

    typedef __attribute__((aligned(64))) struct record {
        int segs;
        int write;
    } record;


    void show(float *c) {
        int i;
        for (i = 0; i < 256; i++)
            if (i % 16 == 0) cout << " [" << i << "] = " << c[i] << endl;
    }

    double getAbs(double x) { return (x < 0) ? (-x) : (x); }

    inline int coordcmp(const void *v1, const void *v2) {
        struct Coordinate *c1 = (struct Coordinate *) v1;
        struct Coordinate *c2 = (struct Coordinate *) v2;

        if (c1->x != c2->x) {
            return (c1->x - c2->x);
        } else {
            return (c1->y - c2->y);
        }
    }

    inline double get_simd(__m512d a, int segs_xx) {
        __mmask8 k_pd, k_n_pd;

        k_pd = _mm512_int2mask(0x01 << segs_xx);
        //    k_n_pd = _mm512_knot(k_pd);

        return _mm512_mask_reduce_add_pd(k_pd, a);
    }

    inline __m512d set_simd_zero(__m512d a, int segs_xx) {
        __mmask8 k_pd, k_n_pd;

        k_pd = _mm512_int2mask(0x01 << segs_xx);
        k_n_pd = _mm512_knot(k_pd);

        __m512d v_zero_pd = _mm512_set1_pd(0);

        return _mm512_mask_add_pd(v_zero_pd, k_n_pd, a, v_zero_pd);
    }


    void show_int(int *c, int a, int b) {
        int i = 0;
        for (i = a; i < (a + 32); i++)
            if (i % 1 == 0)
                std::cout << "i = " << i << ", c[i] = " << c[i] << endl;
    }


    void show(FLOAT_TYPE *c, int a, int b, int stride) {
        int i = 0;
        for (i = a; i < (a + b); i++)
            if (i % stride == 0)
                std::cout << "i = " << i << ", c[i] = " << c[i] << endl;
        /*
   for ( i = b-256; i< b; i++)
      if ( i%16 ==0 )
         std::cout<<"i = "<<i<<", c[i] = "<<c[i]<<endl;
*/
    }

    void show(FLOAT_TYPE *c, FLOAT_TYPE *d, int a, int b, int stride) {
        int i = 0;
        for (i = a; i < (a + b); i++)
            if (i % stride == 0)
                std::cout << " i = " << i << ":  refOut = " << c[i]
                          << ";  refOut_verify = " << d[i] << endl;
        /*
   for ( i = b-256; i< b; i++)
      if ( i%16 ==0 )
         std::cout<<"i = "<<i<<", c[i] = "<<c[i]<<endl;
*/
    }

    void show_list(int *c, int k) {
        int i = 0, j = 0;
        for (i = 0; i < k; i++) { std::cout << "  " << c[i]; }
        std::cout << endl;
    }

    void show_number(int x) {
        std::cout << "*****************************   " << x << endl;
    }

    void show_int(int *c, int k) {
        int i = 0, j = 0;
        for (i = 0; i < k; i += SIMD_LEN) {
            for (j = 0; j < SIMD_LEN; j++) std::cout << "  " << c[i + j];
            std::cout << endl;
        }
    }

    void show_sub(int *ll, int *lr, int k) {
        int i = 0, j = 0;
        for (i = 0; i < k; i++) {
            std::cout << ll[i] << " " << lr[i] << endl;
            ;
        }
    }

    void show_simd_epi32(__m512i art_tmp) {

        int *tmp = (int *) _mm_malloc(sizeof(int) * 16, 64);
        _mm512_store_si512(tmp, art_tmp);
        int i = 0;
        for (i = 0; i < 16; i++) std::cout << "  [" << i << "] = " << tmp[i];
        std::cout << endl;
        _mm_free(tmp);
    }

    void show_simd_ps(__m512 art_tmp) {

        float *tmp = (float *) _mm_malloc(sizeof(float) * 16, 64);
        _mm512_store_ps(tmp, art_tmp);
        int i = 0;
        for (i = 0; i < 16; i++)
            std::cout << "i = " << i << ", tmp[i] = " << tmp[i] << endl;
    }

    void show_simd_pd(__m512d art_tmp) {

        double *tmp = (double *) _mm_malloc(sizeof(double) * 8, 64);
        _mm512_store_pd(tmp, art_tmp);
        int i = 0;
        for (i = 0; i < 8; i++)
            std::cout << "i = " << i << ", tmp[i] = " << tmp[i] << endl;
    }


    // ****************************************************************************
    // Function: readMatrix
    //
    // Purpose:
    //   Reads a sparse matrix from a file of Matrix Market format
    //   Returns the data structures for the CSR format
    //
    // Arguments:
    //   filename: c string with the name of the file to be opened
    //   val_ptr: input - pointer to uninitialized pointer
    //            output - pointer to array holding the non-zero values
    //                     for the  matrix
    //   cols_ptr: input - pointer to uninitialized pointer
    //             output - pointer to array of column indices for each
    //                      element of the sparse matrix
    //   rowDelimiters: input - pointer to uninitialized pointer
    //                  output - pointer to array holding
    //                           indices to rows of the matrix
    //   n: input - pointer to uninitialized int
    //      output - pointer to an int holding the number of non-zero
    //               elements in the matrix
    //   size: input - pointer to uninitialized int
    //         output - pointer to an int holding the number of rows in
    //                  the matrix
    //
    // Programmer: Lukasz Wesolowski
    // Creation: July 2, 2010
    // Returns:  nothing directly
    //           allocates and returns *val_ptr, *cols_ptr, and
    //           *rowDelimiters_ptr indirectly
    //           returns n and size indirectly through pointers
    // ****************************************************************************
    //template <typename floatType>
    void readMatrix(char *filename, floatType **val_ptr, int **cols_ptr,
                    int **rowDelimiters_ptr, int *n, int *numRows,
                    int *numCols) {
        std::string line;
        char id[FIELD_LENGTH];
        char object[FIELD_LENGTH];
        char format[FIELD_LENGTH];
        char field[FIELD_LENGTH];
        char symmetry[FIELD_LENGTH];

        std::ifstream mfs(filename);
        if (!mfs.good()) {
            std::cerr << "Error: unable to open matrix file " << filename
                      << std::endl;
            exit(1);
        }

        int symmetric = 0;
        int pattern = 0;
        int field_complex = 0;

        int nRows, nCols, nElements;

        struct Coordinate *coords;

        // read matrix header
        if (getline(mfs, line).eof()) {
            std::cerr << "Error: file " << filename
                      << " does not store a matrix" << std::endl;
            exit(1);
        }


        sscanf(line.c_str(), "%s %s %s %s %s", id, object, format, field,
               symmetry);

        if (strcmp(object, "matrix") != 0) {
            fprintf(stderr, "Error: file %s does not store a matrix\n",
                    filename);
            exit(1);
        }

        if (strcmp(format, "coordinate") != 0) {
            fprintf(stderr, "Error: matrix representation is dense\n");
            exit(1);
        }

        if (strcmp(field, "pattern") == 0) { pattern = 1; }

        if (strcmp(field, "complex") == 0) { field_complex = 1; }

        if (strcmp(symmetry, "symmetric") == 0) { symmetric = 1; }

        //    pattern = 1;
        //    symmetric = 1;


        while (!getline(mfs, line).eof()) {
            if (line[0] != '%') { break; }
        }

        // read the matrix size and number of non-zero elements
        sscanf(line.c_str(), "%d %d %d", &nRows, &nCols, &nElements);
        //    sscanf(line.c_str(), "%d %d", &nRows, &nElements);
        //    nCols = nRows;

        int nElements_padding =
                (nElements % 16 == 0) ? nElements : (nElements + 16) / 16 * 16;

        //    std::cout<<" numRows is "<<nRows<<";  Number of Elements is "<< nElements<<";  After Padding, Number of Elements is "<< nElements_padding<<endl;

        int valSize = nElements_padding * sizeof(struct Coordinate);

        if (symmetric) { valSize *= 2; }

        //    coords = new Coordinate[valSize];
        coords = (struct Coordinate *) malloc(valSize);


        //    std::cout<<"11111111111111111111111111111"<<endl;


        int index = 0;
        float xx99 = 0;
        while (!getline(mfs, line).eof()) {
            if (pattern) {
                sscanf(line.c_str(), "%d %d", &coords[index].x,
                       &coords[index].y);
                //            coords[index].val = 1;
                coords[index].val = index % 13;
                //            coords[index].y = 1;
                // assign a random value
                //            coords[index].val = ((floatType) MAX_RANDOM_VAL *
                //                                 (rand() / (RAND_MAX + 1.0)));
            } else if (field_complex) {
                // read the value from file
                sscanf(line.c_str(), "%d %d %f %f", &coords[index].x,
                       &coords[index].y, &coords[index].val, xx99);
            } else {
                // read the value from file
                sscanf(line.c_str(), "%d %d %f", &coords[index].x,
                       &coords[index].y, &coords[index].val);
            }

            // convert into index-0-as-start representation
            //        coords[index].x--;
            //        coords[index].y--;
            index++;
            //        std::cout<<"index = "<<index<<endl;

            // add the mirror element if not on main diagonal
            if (symmetric && coords[index - 1].x != coords[index - 1].y) {
                coords[index].x = coords[index - 1].y;
                coords[index].y = coords[index - 1].x;
                coords[index].val = coords[index - 1].val;
                index++;
            }
        }

        //    std::cout<<"22222222222222222222222222222222222"<<endl;

        nElements = index;

        nElements_padding =
                (nElements % 16 == 0) ? nElements : (nElements + 16) / 16 * 16;

        std::cout << "========================================================="
                     "=================="
                  << endl;
        std::cout << "=========*********  Informations of the sparse matrix   "
                     "*********=========="
                  << endl;
        std::cout << endl;
        //    std::cout<<" numRows is "<<nRows<<"; numCols is "<< nCols<<";  Number of Elements is "<< nElements<<";  After Padding, Number of Elements is "<< nElements_padding<<endl;
        std::cout << "     Number of Rows is :" << nRows << endl;
        std::cout << "  Number of Columns is :" << nCols << endl;
        std::cout << " Number of Elements is :" << nElements << endl;
        std::cout << "       After Alignment :" << nElements_padding << endl;
        std::cout << endl;
        std::cout << "========================================================="
                     "=================="
                  << endl;

        std::cout << "............ Converting the Raw matrix to CSR "
                     "................."
                  << endl;

        //    std::cout<<" index = "<< index <<";  nElements_padding = "<< nElements_padding<<endl;

        for (int qq = index; qq < nElements_padding; qq++) {
            coords[qq].x = coords[index - 1].x;
            coords[qq].y = coords[index - 1].y;
            coords[qq].val = 0;

            //            std::cout<<"Padding: qq = "<< qq<<";  x = "<< coords[qq].x <<";  y = "<< coords[qq].y <<"; val = "<< coords[qq].val <<endl;
        }

        //sort the elements
        qsort(coords, nElements_padding, sizeof(struct Coordinate), coordcmp);

        //    std::cout<<"33333333333333333333333333333333333"<<endl;

        // create CSR data structures
        *n = nElements_padding;
        *numRows = nRows;
        *numCols = nCols;
        //    *val_ptr = new floatType[nElements_padding];
        //    *cols_ptr = new int[nElements_padding];
        //    *rowDelimiters_ptr = new int[nRows+2];

        *val_ptr = (floatType *) _mm_malloc(
                sizeof(floatType) * nElements_padding, 64);
        *cols_ptr = (int *) _mm_malloc(sizeof(int) * nElements_padding, 64);
        *rowDelimiters_ptr = (int *) _mm_malloc(sizeof(int) * (nRows + 2), 64);

        floatType *val = *val_ptr;
        int *cols = *cols_ptr;
        int *rowDelimiters = *rowDelimiters_ptr;

        rowDelimiters[0] = 0;
        //    rowDelimiters[nRows] = nElements_padding;
        int r = 0;

        //    std::cout<<"444444444444444444444444444444444444"<<endl;
        int i = 0;
        for (i = 0; i < nElements_padding; i++) {
            while (coords[i].x != r) { rowDelimiters[++r] = i; }
            val[i] = coords[i].val;
            cols[i] = coords[i].y;
            //        std::cout<<"i = "<<i<<endl;
        }

        for (int k = r + 1; k <= (nRows + 1); k++) {
            rowDelimiters[k] = i - 1;
            //       std::cout<<" rowDelimiter["<<k<<"] = "<< rowDelimiters[k]<<endl;
        }

        r = 0;

        //    delete[] coords;

        free(coords);
    }


    // ****************************************************************************
    // Function: fill
    //
    // Purpose:
    //   Simple routine to initialize input array
    //
    // Arguments:
    //   A: pointer to the array to initialize
    //   n: number of elements in the array
    //   maxi: specifies range of random values
    //
    // Programmer: Lukasz Wesolowski
    // Creation: June 21, 2010
    // Returns:  nothing
    //
    // ****************************************************************************
    //template <typename floatType>
    void fill(floatType *A, const int n, const float maxi) {
        for (int j = 0; j < n; j++) { A[j] = 1.0; }
    }

    void pre_processing(int Nthrds, int N_start, int N_step, int *vPack_Nblock,
                        int *vPack_vec_record, int *vPack_nnz_rows,
                        floatType *vPack_vec_vals, int *vPack_vec_cols,
                        floatType *h_val, int *h_cols, int *vPack_vec_final,
                        int *vPack_vec_final_2, int *vPack_split,
                        floatType *refOut, int nItems, int numRows, int omega,
                        int *h_rowDelimiters, char *filename) {
#pragma omp parallel num_threads(Nthrds)
        {
            int thread_idx = omp_get_thread_num();

            int thread_nnz = (nItems / Nthrds / omega / 16) * 16 * omega;

            int thread_break = (nItems - thread_nnz * Nthrds) / 16;

            int thread_nnz_2 = (nItems / Nthrds / omega / 16 + 1) * 16;

            if (thread_idx >= N_start && thread_idx < (N_start + N_step)) {
                int seg_idx = 0;

                for (int row_block = 0; row_block < vPack_Nblock[thread_idx];
                     row_block++) {

                    if (seg_idx != 0) seg_idx = (seg_idx / 16 + 1) * 16;


                    int thread_nnz_start;
                    int thread_nnz_end;

                    if (thread_idx < thread_break) {
                        thread_nnz_start = thread_idx * (thread_nnz + 16);
                        thread_nnz_end = (thread_idx + 1) * (thread_nnz + 16);
                    } else {
                        thread_nnz_start =
                                thread_idx * thread_nnz + thread_break * 16;
                        thread_nnz_end = (thread_idx + 1) * thread_nnz +
                                         thread_break * 16;
                    }

                    if (thread_idx == (Nthrds - 1) &&
                        (row_block == (vPack_Nblock[thread_idx] - 1) ||
                         vPack_Nblock[thread_idx] == 1))
                        thread_nnz_end = nItems;

                    int thread_rows_start, thread_rows_end;

                    int start = 0;
                    int stop = numRows;
                    int median;
                    int key_median;
                    int key_median_2;

                    while (stop >= start) {
                        median = (stop + start) / 2;

                        key_median = h_rowDelimiters[median];
                        key_median_2 = h_rowDelimiters[median + 1];

                        if (thread_nnz_start >= key_median) start = median + 1;
                        else
                            stop = median - 1;
                    }

                    thread_rows_start = start - 1;

                    start = thread_rows_start;
                    stop = numRows;

                    while (stop >= start) {

                        median = (stop + start) / 2;
                        key_median = h_rowDelimiters[median];

                        if ((thread_nnz_end - 1) >= key_median)
                            start = median + 1;
                        else
                            stop = median - 1;
                    }

                    thread_rows_end = start - 1;

                    while ((h_rowDelimiters[thread_rows_end + 1] -
                                    h_rowDelimiters[thread_rows_end] ==
                            0) &&
                           (thread_rows_end <= numRows))
                        thread_rows_end++;

                    vPack_nnz_rows[thread_idx * 4] = thread_nnz_start;
                    vPack_nnz_rows[thread_idx * 4 + 1] = thread_nnz_end;

                    vPack_nnz_rows[thread_idx * 4 + 2] = thread_rows_start;
                    vPack_nnz_rows[thread_idx * 4 + 3] = thread_rows_end;

                    int thread_rows_span =
                            thread_rows_end - thread_rows_start + 1;

                    floatType *vPack_vals = vPack_vec_vals + thread_nnz_start;
                    int *vPack_cols = vPack_vec_cols + thread_nnz_start;

                    floatType *orig_vals = h_val + thread_nnz_start;
                    int *orig_cols = h_cols + thread_nnz_start;

                    int *vPack_final = vPack_vec_final + thread_idx * 16;
                    int *vPack_final_2 = vPack_vec_final_2 +
                                         thread_idx * 16 * omega +
                                         16 * row_block;
                    int *vPack_record =
                            vPack_vec_record +
                            2 * (thread_idx * 32 + thread_rows_start) / 16 * 16;

                    zmmd_t vPack_valID, vPack_rowID, vPack_count;
                    vPack_count.regi32 = _mm512_set1_epi32(1);
                    vPack_valID.regi32 = _mm512_set1_epi32(0);
                    vPack_rowID.regi32 = _mm512_set1_epi32(0);

                    zmmd_t vPack_flag;
                    vPack_flag.regi32 = _mm512_set1_epi32(-1);

                    int thread_rows_start_init = thread_rows_start;


                    for (int i = 0; i < 8; i++) {
                        if (thread_rows_start < thread_rows_end) {
                            vPack_valID.elemsi32[i] =
                                    h_rowDelimiters[thread_rows_start] -
                                    thread_nnz_start;
                            vPack_rowID.elemsi32[i] = thread_rows_start;
                            vPack_count.elemsi32[i] =
                                    h_rowDelimiters[thread_rows_start + 1] -
                                    h_rowDelimiters[thread_rows_start];
                        } else if (thread_rows_start == thread_rows_end) {
                            vPack_valID.elemsi32[i] =
                                    h_rowDelimiters[thread_rows_start] -
                                    thread_nnz_start;
                            vPack_rowID.elemsi32[i] = thread_rows_start;
                            vPack_count.elemsi32[i] =
                                    thread_nnz_end -
                                    h_rowDelimiters[thread_rows_start];
                        } else if (thread_rows_start > thread_rows_end) {
                            vPack_valID.elemsi32[i] = 0;
                            vPack_rowID.elemsi32[i] = 0;
                            vPack_count.elemsi32[i] = 0;
                        }

                        if (i == 0) {
                            vPack_valID.elemsi32[i] = 0;
                            vPack_count.elemsi32[i] =
                                    h_rowDelimiters[thread_rows_start + 1] -
                                    thread_nnz_start;
                            if (thread_rows_start == thread_rows_end) {
                                vPack_count.elemsi32[i] =
                                        thread_nnz_end - thread_nnz_start;
                            }
                        }

                        thread_rows_start++;//should after the if(i==0)
                    }

                    __m512i v_zero_epi32 = _mm512_set1_epi32(0);
                    __m512i v_one_epi32 = _mm512_set1_epi32(1);
                    __m512i v_two_epi32 = _mm512_set1_epi32(2);
                    __m512i v_four_epi32 = _mm512_set1_epi32(4);
                    __m512i v_eight_epi32 = _mm512_set1_epi32(8);
                    __m512i v_sixteen_epi32 = _mm512_set1_epi32(16);
                    __m512i v_0o8_epi32 =
                            _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6,
                                             5, 4, 3, 2, 1, 0);
                    __m512d v_zero_pd = _mm512_set1_pd(0);
                    __m512d v_one_pd = _mm512_set1_pd(1);

                    __m512d v_vals = _mm512_set1_pd(0);
                    __m512i v_cols, v_trans;

                    __mmask16 k_lo = _mm512_int2mask(0x00ff);
                    __mmask16 k_hi = _mm512_int2mask(0xff00);

                    __mmask16 mask_0xAAAA = _mm512_int2mask(0xaaaa);
                    __mmask16 mask_0x5555 = _mm512_int2mask(0x5555);


                    int i = 0;

                    int first_in = 0;
                    int tt_index = 0;

                    zmmd_t tt_record;
                    tt_record.regi32 = _mm512_set1_epi32(0);

                    __m512i v_idx_cols =
                            _mm512_set_epi32(15, 7, 13, 5, 11, 3, 9, 1, 14, 6,
                                             12, 4, 10, 2, 8, 0);

                    zmmd_t first_flag;
                    first_flag.regi32 = _mm512_set1_epi32(0);

                    int kk = 0;

                    int kkdu = 0;

#pragma vector nontemporal(vPack_vals, vPack_cols, vPack_record)
                    for (i = 0; i < (thread_nnz_end - thread_nnz_start) / 8;
                         i++) {
                        int xx = _mm512_mask_reduce_min_epi32(
                                k_lo, vPack_count.regi32);

                        if (xx == 0) {
                            for (int kk = 0; kk < 8; kk++) {
                                if (vPack_count.elemsi32[kk] == 0) {
                                    if (thread_rows_start <= thread_rows_end) {
                                        if (vPack_rowID.elemsi32[kk] ==
                                            thread_rows_start_init) {
                                            vPack_split[thread_idx * 2] =
                                                    i * 8 + kk;
                                        } else {
                                            vPack_record[seg_idx] = i * 8 + kk;
                                            vPack_record[seg_idx + 1] =
                                                    vPack_rowID.elemsi32[kk];
                                            seg_idx += 2;
                                        }

                                        while (h_rowDelimiters[thread_rows_start +
                                                               1] -
                                                       h_rowDelimiters
                                                               [thread_rows_start] ==
                                               0)
                                            thread_rows_start++;

                                        vPack_valID.elemsi32[kk] =
                                                h_rowDelimiters
                                                        [thread_rows_start] -
                                                thread_nnz_start;
                                        vPack_rowID.elemsi32[kk] =
                                                thread_rows_start;
                                        vPack_count.elemsi32[kk] =
                                                h_rowDelimiters
                                                        [thread_rows_start +
                                                         1] -
                                                h_rowDelimiters
                                                        [thread_rows_start];

                                        if (thread_rows_start ==
                                            thread_rows_end) {
                                            if (vPack_split[thread_idx * 2 +
                                                            1] == 0) {
                                                vPack_split[thread_idx * 2 +
                                                            1] = i * 8 + kk;
                                            }

                                            vPack_count.elemsi32[kk] =
                                                    thread_nnz_end -
                                                    h_rowDelimiters
                                                            [thread_rows_start];

                                            _mm512_mask_store_epi32(
                                                    vPack_final_2, k_lo,
                                                    vPack_rowID.regi32);

                                            __mmask16 kzz =
                                                    _mm512_cmpeq_epi32_mask(
                                                            vPack_count.regi32,
                                                            v_zero_epi32);
                                            vPack_flag.regi32 =
                                                    _mm512_mask_add_epi32(
                                                            vPack_flag.regi32,
                                                            kzz, v_zero_epi32,
                                                            v_zero_epi32);
                                        }
                                        thread_rows_start++;
                                    } else if (thread_rows_start >
                                               thread_rows_end) {
                                        int ave = _mm512_mask_reduce_add_epi32(
                                                          k_lo,
                                                          vPack_count.regi32) /
                                                  8;
                                        int cali_i = 0;
                                        for (cali_i = 0; cali_i < 8; cali_i++)
                                            if (vPack_count.elemsi32[cali_i] >
                                                ave)
                                                break;

                                        if (first_flag.elemsi32[kk] == 0) {
                                            if (first_in == 0) {
                                                if (vPack_split[thread_idx * 2 +
                                                                1] == 0) {
                                                    if (thread_rows_span <= 8)
                                                        vPack_split[thread_idx *
                                                                            2 +
                                                                    1] = -1;
                                                    else
                                                        vPack_split[thread_idx *
                                                                            2 +
                                                                    1] =
                                                                i * 8 + kk;
                                                }

                                                _mm512_mask_store_epi32(
                                                        vPack_final_2, k_lo,
                                                        vPack_rowID.regi32);

                                                first_in = 1;
                                            }

                                            vPack_record[seg_idx] = i * 8 + kk;
                                            vPack_record[seg_idx + 1] = kk;
                                            vPack_flag.elemsi32[kk] = cali_i;

                                            first_flag.elemsi32[kk] = 1;
                                        } else {
                                            vPack_record[seg_idx] = i * 8 + kk;
                                            vPack_record[seg_idx + 1] =
                                                    vPack_flag.elemsi32[kk];
                                            vPack_flag.elemsi32[kk] = cali_i;
                                        }
                                        {
                                            vPack_valID.elemsi32[kk] =
                                                    vPack_valID
                                                            .elemsi32[cali_i];
                                            vPack_rowID.elemsi32[kk] = cali_i;
                                            vPack_count.elemsi32[kk] = ave;
                                            vPack_count.elemsi32[cali_i] =
                                                    vPack_count
                                                            .elemsi32[cali_i] -
                                                    ave;
                                            vPack_valID.elemsi32[cali_i] =
                                                    vPack_valID
                                                            .elemsi32[cali_i] +
                                                    ave;
                                        }
                                        seg_idx += 2;
                                    }
                                }
                            }
                        }
                        v_vals = _mm512_i32logather_pd(vPack_valID.regi32,
                                                       orig_vals, _MM_SCALE_8);
                        _mm512_store_pd(vPack_vals + i * 8, v_vals);

                        if (i % 2 == 0) {
                            v_trans = vPack_valID.regi32;
                        } else {
                            v_trans = _mm512_add_epi32(
                                    _mm512_permute4f128_epi32(
                                            vPack_valID.regi32, _MM_PERM_BADC),
                                    v_trans);
                            v_cols = _mm512_i32gather_epi32(v_trans, orig_cols,
                                                            _MM_SCALE_4);
                            _mm512_store_epi32(vPack_cols + (i - 1) * 8,
                                               v_cols);
                        }

                        vPack_valID.regi32 = _mm512_mask_add_epi32(
                                v_zero_epi32, k_lo, vPack_valID.regi32,
                                v_one_epi32);
                        vPack_count.regi32 = _mm512_mask_sub_epi32(
                                v_zero_epi32, k_lo, vPack_count.regi32,
                                v_one_epi32);

                        if (i == (thread_nnz_end - thread_nnz_start) / 8 - 1) {
                            for (int kk = 0; kk < 8; kk++) {
                                if (vPack_flag.elemsi32[kk] == -1) {
                                    vPack_record[seg_idx] = -1;
                                    vPack_record[seg_idx + 1] = kk;
                                    seg_idx += 2;
                                } else {
                                    vPack_record[seg_idx] = -1;
                                    vPack_record[seg_idx + 1] =
                                            vPack_flag.elemsi32[kk];
                                    seg_idx += 2;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void spmv_compute_kernel(int Nthrds, int N_start, int N_step,
                             int *vPack_Nblock, int *vPack_vec_record,
                             int *vPack_nnz_rows, floatType *vPack_vec_vals,
                             int *vPack_vec_cols, floatType *h_val, int *h_cols,
                             int *vPack_vec_final, int *vPack_vec_final_2,
                             int *vPack_split, floatType *refOut, int nItems,
                             int numRows, int omega, int *h_rowDelimiters,
                             char *filename, floatType *h_vec, int Ntimes) {
#pragma omp parallel num_threads(Nthrds)
        {
            int thread_idx = omp_get_thread_num();


            if (thread_idx >= N_start && thread_idx < (N_start + N_step)) {

                int thread_nnz_2 = (nItems / Nthrds / omega / 16 + 1) * 16;

                int simd_seg_index = 0;

                for (int row_block = 0; row_block < vPack_Nblock[thread_idx];
                     row_block++) {

                    if (simd_seg_index != 0)
                        simd_seg_index = (simd_seg_index / 16 + 1) * 16;


                    int thread_nnz_start = vPack_nnz_rows[thread_idx * 4];
                    int thread_nnz_end = vPack_nnz_rows[thread_idx * 4 + 1];

                    int thread_rows_start = vPack_nnz_rows[thread_idx * 4 + 2];
                    int thread_rows_end = vPack_nnz_rows[thread_idx * 4 + 3];

                    floatType *csr_vals = vPack_vec_vals + thread_nnz_start;
                    int *csr_cols = vPack_vec_cols + thread_nnz_start;

                    int *csr_final = vPack_vec_final + thread_idx * 16;
                    int *csr_final_2 = vPack_vec_final_2 +
                                       thread_idx * 16 * omega + row_block * 16;
                    int *csr_record =
                            vPack_vec_record +
                            2 * (thread_idx * 32 + thread_rows_start) / 16 * 16;

                    zmmd_t z_rets;
                    z_rets.reg = _mm512_set1_pd(0);

                    __m512d r_rets;
                    r_rets = _mm512_set1_pd(0);

                    zmmd_t t_rets;
                    t_rets.reg = _mm512_set1_pd(0);

                    int t_index = 0;

                    __m512i v_zero_epi32 = _mm512_set1_epi32(0);
                    __m512d v_zero_pd = _mm512_set1_pd(0);
                    __m512d v_one_pd = _mm512_set1_pd(1);

                    __m512i v_two_epi32 = _mm512_set1_epi32(2);
                    __m512i v_single_epi32 = _mm512_set_epi32(
                            0, 0, 0, 0, 0, 0, 0, 0, 15, 13, 11, 9, 7, 5, 3, 1);

                    __mmask16 k_lo = _mm512_int2mask(0x00ff);
                    __mmask16 k_hi = _mm512_int2mask(0xff00);

                    __m512d v_vec, v_vals;
                    __m512i v_cols;


                    int segs_xx;

                    int ncsr_start = vPack_split[thread_idx * 2];
                    int ncsr = vPack_split[thread_idx * 2 + 1];

                    __m512i v_write;

                    int prefetch_distance = 8;

                    int i = 0;

                    int ii = 0;

#pragma noprefetch h_vec
                    for (i = 0; i < ncsr_start / 8 * 8;
                         i += prefetch_distance) {
                        int *prefetch = &csr_cols[i + prefetch_distance];
                        for (int k = 0; k < prefetch_distance / 8; k++) {
                            int idx0 = prefetch[k * 8 + 0];
                            int idx1 = prefetch[k * 8 + 1];
                            int idx2 = prefetch[k * 8 + 2];
                            int idx3 = prefetch[k * 8 + 3];
                            int idx4 = prefetch[k * 8 + 4];
                            int idx5 = prefetch[k * 8 + 5];
                            int idx6 = prefetch[k * 8 + 6];
                            int idx7 = prefetch[k * 8 + 7];


                            _mm_prefetch((const char *) &h_vec[idx0],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx1],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx2],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx3],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx4],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx5],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx6],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx7],
                                         _MM_HINT_T1);
                        }

                        for (int ii = 0; ii < prefetch_distance; ii += 8) {

                            while ((csr_record[simd_seg_index] & 0xfffffff8) ==
                                   ((i + ii) & 0xfffffff8)) {
                                segs_xx =
                                        csr_record[simd_seg_index] & 0x00000007;

                                int row_idd = csr_record[simd_seg_index + 1];
                                refOut[row_idd] = get_simd(r_rets, segs_xx);
                                r_rets = set_simd_zero(r_rets, segs_xx);

                                simd_seg_index += 2;
                            }
                            v_cols = ((i + ii) & 0x0000000f) == 0
                                             ? _mm512_load_epi32(csr_cols +
                                                                 (i + ii))
                                             : _mm512_permute4f128_epi32(
                                                       v_cols, _MM_PERM_BADC);
                            v_vec = _mm512_i32logather_pd(v_cols, h_vec,
                                                          _MM_SCALE_8);
                            v_vals = _mm512_load_pd(csr_vals + i + ii);
                            r_rets = _mm512_fmadd_pd(v_vals, v_vec, r_rets);
                        }
                    }


#pragma noprefetch h_vec
                    for (i = ncsr_start / 8 * 8; i < ncsr_start / 8 * 8 + 8;
                         i += prefetch_distance) {
                        int *prefetch = &csr_cols[i + prefetch_distance];
                        for (int k = 0; k < prefetch_distance / 8; k++) {
                            int idx0 = prefetch[k * 8 + 0];
                            int idx1 = prefetch[k * 8 + 1];
                            int idx2 = prefetch[k * 8 + 2];
                            int idx3 = prefetch[k * 8 + 3];
                            int idx4 = prefetch[k * 8 + 4];
                            int idx5 = prefetch[k * 8 + 5];
                            int idx6 = prefetch[k * 8 + 6];
                            int idx7 = prefetch[k * 8 + 7];


                            _mm_prefetch((const char *) &h_vec[idx0],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx1],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx2],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx3],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx4],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx5],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx6],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx7],
                                         _MM_HINT_T1);
                        }

                        // because the link is padded with cols 0; vals 0 at last. there may came across some errors with elements[0]
                        for (int ii = 0; ii < prefetch_distance; ii += 8) {

#pragma omp atomic
                            refOut[thread_rows_start] +=
                                    get_simd(r_rets, ncsr_start % 8);
                            r_rets = set_simd_zero(r_rets, ncsr_start % 8);


                            while ((csr_record[simd_seg_index] & 0xfffffff8) ==
                                   ((i + ii) & 0xfffffff8)) {
                                segs_xx =
                                        csr_record[simd_seg_index] & 0x00000007;


                                int row_idd = csr_record[simd_seg_index + 1];
                                refOut[row_idd] = get_simd(r_rets, segs_xx);
                                r_rets = set_simd_zero(r_rets, segs_xx);
                                simd_seg_index += 2;
                            }

                            v_cols = ((i + ii) & 0x0000000f) == 0
                                             ? _mm512_load_epi32(csr_cols +
                                                                 (i + ii))
                                             : _mm512_permute4f128_epi32(
                                                       v_cols, _MM_PERM_BADC);
                            v_vec = _mm512_i32logather_pd(v_cols, h_vec,
                                                          _MM_SCALE_8);
                            v_vals = _mm512_load_pd(csr_vals + i + ii);
                            r_rets = _mm512_fmadd_pd(v_vals, v_vec, r_rets);
                        }
                    }

                    z_rets.reg = r_rets;


#pragma noprefetch h_vec
                    for (i = ncsr_start / 8 * 8 + 8; i < ncsr / 8 * 8;
                         i += prefetch_distance) {
                        int *prefetch = &csr_cols[i + prefetch_distance];
                        for (int k = 0; k < prefetch_distance / 8; k++) {
                            int idx0 = prefetch[k * 8 + 0];
                            int idx1 = prefetch[k * 8 + 1];
                            int idx2 = prefetch[k * 8 + 2];
                            int idx3 = prefetch[k * 8 + 3];
                            int idx4 = prefetch[k * 8 + 4];
                            int idx5 = prefetch[k * 8 + 5];
                            int idx6 = prefetch[k * 8 + 6];
                            int idx7 = prefetch[k * 8 + 7];


                            _mm_prefetch((const char *) &h_vec[idx0],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx1],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx2],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx3],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx4],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx5],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx6],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx7],
                                         _MM_HINT_T1);
                        }

                        for (int ii = 0; ii < prefetch_distance; ii += 8) {
                            while ((csr_record[simd_seg_index] & 0xfffffff8) ==
                                   ((i + ii) & 0xfffffff8)) {
                                segs_xx =
                                        csr_record[simd_seg_index] & 0x00000007;

                                int row_idd = csr_record[simd_seg_index + 1];
                                refOut[row_idd] = z_rets.elems[segs_xx];
                                z_rets.elems[segs_xx] = 0;

                                simd_seg_index += 2;
                            }

                            v_cols = ((i + ii) & 0x0000000f) == 0
                                             ? _mm512_load_epi32(csr_cols +
                                                                 (i + ii))
                                             : _mm512_permute4f128_epi32(
                                                       v_cols, _MM_PERM_BADC);

                            v_vec = _mm512_i32logather_pd(v_cols, h_vec,
                                                          _MM_SCALE_8);
                            v_vals = _mm512_load_pd(csr_vals + i + ii);
                            z_rets.reg =
                                    _mm512_fmadd_pd(v_vals, v_vec, z_rets.reg);
                        }
                    }
                    if (ncsr >= 0)
#pragma noprefetch h_vec
                        for (i = ncsr / 8 * 8; i < ncsr / 8 * 8 + 8;
                             i += prefetch_distance) {
                            int *prefetch = &csr_cols[i + prefetch_distance];
                            for (int k = 0; k < prefetch_distance / 8; k++) {
                                int idx0 = prefetch[k * 8 + 0];
                                int idx1 = prefetch[k * 8 + 1];
                                int idx2 = prefetch[k * 8 + 2];
                                int idx3 = prefetch[k * 8 + 3];
                                int idx4 = prefetch[k * 8 + 4];
                                int idx5 = prefetch[k * 8 + 5];
                                int idx6 = prefetch[k * 8 + 6];
                                int idx7 = prefetch[k * 8 + 7];


                                _mm_prefetch((const char *) &h_vec[idx0],
                                             _MM_HINT_T1);
                                _mm_prefetch((const char *) &h_vec[idx1],
                                             _MM_HINT_T1);
                                _mm_prefetch((const char *) &h_vec[idx2],
                                             _MM_HINT_T1);
                                _mm_prefetch((const char *) &h_vec[idx3],
                                             _MM_HINT_T1);
                                _mm_prefetch((const char *) &h_vec[idx4],
                                             _MM_HINT_T1);
                                _mm_prefetch((const char *) &h_vec[idx5],
                                             _MM_HINT_T1);
                                _mm_prefetch((const char *) &h_vec[idx6],
                                             _MM_HINT_T1);
                                _mm_prefetch((const char *) &h_vec[idx7],
                                             _MM_HINT_T1);
                            }

                            for (int ii = 0; ii < prefetch_distance; ii += 8) {
                                while ((csr_record[simd_seg_index] &
                                        0xfffffff8) ==
                                       ((i + ii) & 0xfffffff8)) {
                                    segs_xx = csr_record[simd_seg_index] &
                                              0x00000007;
                                    if (csr_record[simd_seg_index] <= ncsr) {
                                        int row_idd =
                                                csr_record[simd_seg_index + 1];
                                        refOut[row_idd] = z_rets.elems[segs_xx];
                                        z_rets.elems[segs_xx] = 0;

                                        simd_seg_index += 2;
                                    } else {
                                        if (t_index != -1) {
                                            t_index = -1;
                                            t_rets.reg = v_zero_pd;

                                            t_rets.elems[segs_xx] +=
                                                    z_rets.elems[segs_xx];
                                            z_rets.elems[segs_xx] = 0;
                                        } else {
                                            int tmp_idx =
                                                    csr_record[simd_seg_index +
                                                               1];
                                            t_rets.elems[tmp_idx] +=
                                                    z_rets.elems[segs_xx];
                                            z_rets.elems[segs_xx] = 0;
                                            simd_seg_index += 2;
                                        }
                                    }
                                }

                                if (t_index >= 0)
                                    for (int iik = 0; iik < t_index; iik++) {
                                        int tmp_w = csr_record[simd_seg_index -
                                                               t_index * 2 +
                                                               iik * 2 + 1];
                                        refOut[tmp_w] += t_rets.elems[iik];
                                        t_rets.elems[iik] = 0;
                                    }


                                v_cols = ((i + ii) & 0x0000000f) == 0
                                                 ? _mm512_load_epi32(csr_cols +
                                                                     (i + ii))
                                                 : _mm512_permute4f128_epi32(
                                                           v_cols,
                                                           _MM_PERM_BADC);

                                v_vec = _mm512_i32logather_pd(v_cols, h_vec,
                                                              _MM_SCALE_8);
                                v_vals = _mm512_load_pd(csr_vals + i + ii);
                                z_rets.reg = _mm512_fmadd_pd(v_vals, v_vec,
                                                             z_rets.reg);
                            }
                        }


#pragma noprefetch h_vec
                    for (i = ncsr / 8 * 8 + 8;
                         i < (thread_nnz_end - thread_nnz_start);
                         i += prefetch_distance) {
                        int *prefetch = &csr_cols[i + prefetch_distance];
                        for (int k = 0; k < prefetch_distance / 8; k++) {
                            int idx0 = prefetch[k * 8 + 0];
                            int idx1 = prefetch[k * 8 + 1];
                            int idx2 = prefetch[k * 8 + 2];
                            int idx3 = prefetch[k * 8 + 3];
                            int idx4 = prefetch[k * 8 + 4];
                            int idx5 = prefetch[k * 8 + 5];
                            int idx6 = prefetch[k * 8 + 6];
                            int idx7 = prefetch[k * 8 + 7];


                            _mm_prefetch((const char *) &h_vec[idx0],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx1],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx2],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx3],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx4],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx5],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx6],
                                         _MM_HINT_T1);
                            _mm_prefetch((const char *) &h_vec[idx7],
                                         _MM_HINT_T1);
                        }

                        // because the link is padded with cols 0; vals 0 at last. there may came across some errors with elements[0]
                        for (int ii = 0; ii < prefetch_distance; ii += 8) {
                            while ((csr_record[simd_seg_index] & 0xfffffff8) ==
                                   ((i + ii) & 0xfffffff8)) {
                                segs_xx =
                                        csr_record[simd_seg_index] & 0x00000007;

                                int tmp_idx = csr_record[simd_seg_index + 1];

                                t_rets.elems[tmp_idx] += z_rets.elems[segs_xx];
                                z_rets.elems[segs_xx] = 0;
                                simd_seg_index += 2;
                            }

                            v_cols = ((i + ii) & 0x0000000f) == 0
                                             ? _mm512_load_epi32(csr_cols +
                                                                 (i + ii))
                                             : _mm512_permute4f128_epi32(
                                                       v_cols, _MM_PERM_BADC);

                            v_vec = _mm512_i32logather_pd(v_cols, h_vec,
                                                          _MM_SCALE_8);
                            v_vals = _mm512_load_pd(csr_vals + i + ii);
                            z_rets.reg =
                                    _mm512_fmadd_pd(v_vals, v_vec, z_rets.reg);
                        }
                    }

                    if (i == (thread_nnz_end - thread_nnz_start)) {
                        for (int ttk = 0; ttk < 8; ttk++) {
                            int tmp_idx = csr_record[simd_seg_index + 1];
                            simd_seg_index += 2;
                            t_rets.elems[tmp_idx] += z_rets.elems[ttk];
                        }

                        for (int ttk = 0; ttk < 8; ttk++) {
                            int tmp_w = csr_final_2[ttk];
#pragma omp atomic
                            refOut[tmp_w] += t_rets.elems[ttk];
                        }
                    }
                }
            }
        }
    }

    void build_set(ExecutorParameters& ep) {
#ifdef MMAP
        h_vec = (floatType *) mmap(
                0, numCols * sizeof(double), PROT_READ | PROT_WRITE,
                MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB, -1, 0);
        //    hbw_posix_memalign_psize(&h_vec, 64, sizeof(floatType) * numCols, HBW_PAGESIZE_2MB);
#else
#endif

        struct v_Range *range_thread_nnz =
                (struct v_Range *) malloc(sizeof(struct v_Range) * ep.Nthrds);

        ep.vPack_vec_vals =
                (floatType *) _mm_malloc(sizeof(floatType) * ep.nItems, 64);
        ep.vPack_vec_cols = (int *) _mm_malloc(sizeof(int) * ep.nItems, 64);

        ep.vPack_vec_final =
                (int *) _mm_malloc(sizeof(int) * 16 * ep.Nthrds, 64);
        ep.vPack_vec_final_2 =
                (int *) _mm_malloc(sizeof(int) * 16 * ep.omega * ep.Nthrds, 64);

        ep.vPack_vec_record = (int *) _mm_malloc(
                sizeof(int) * 2 * (ep.numRows + 240 + ep.Nthrds * 32), 64);

        if (ep.vPack_vec_record == NULL)
            std::cout << " vPack_vec_record is NULL "
                         "***********************************"
                      << endl;


        ep.vPack_split = (int *) malloc(sizeof(int) * ep.Nthrds *
                                          2);// must be inited to be 0

                                          ep.vPack_Nblock =
                (int *) malloc(sizeof(int) * ep.Nthrds);// must be inited to be 0
                ep.vPack_nnz_rows = (int *) _mm_malloc(sizeof(int) * 4 * ep.Nthrds,
                                                 64);// must be inited to be 0

                                                 for (int i = 0; i < ep.Nthrds * 2; i++) ep.vPack_split[i] = 0;
                                                 for (int i = 0; i < ep.Nthrds; i++)
                                                     if (i < ep.Nthrds) ep.vPack_Nblock[i] = 1;
            else
                ep.vPack_Nblock[i] = ep.omega;

#pragma omp parallel num_threads(ep.Nthrds)
for (int i = 0; i < ep.numRows; i++) ep.refOut[i] = 0;
            pre_processing(ep.Nthrds, ep.N_start, ep.N_step, ep.vPack_Nblock, ep.vPack_vec_record,
                           ep.vPack_nnz_rows, ep.vPack_vec_vals, ep.vPack_vec_cols, ep.h_val,
                           ep.h_cols, ep.vPack_vec_final, ep.vPack_vec_final_2, ep.vPack_split,
                           ep.refOut, ep.nItems, ep.numRows, ep.omega, ep.h_rowDelimiters,
                           ep.filename);
    }

    void run_spmv(int Nthrds, int N_start, int N_step, int *vPack_Nblock,
                  int *vPack_vec_record, int *vPack_nnz_rows,
                  floatType *vPack_vec_vals, int *vPack_vec_cols,
                  floatType *h_val, int *h_cols, int *vPack_vec_final,
                  int *vPack_vec_final_2, int *vPack_split, floatType *refOut,
                  int nItems, int numRows, int omega, int *h_rowDelimiters,
                  char *filename, floatType *h_vec, int Ntimes) {
        spmv_compute_kernel(
                Nthrds, N_start, N_step, vPack_Nblock, vPack_vec_record,
                vPack_nnz_rows, vPack_vec_vals, vPack_vec_cols, h_val, h_cols,
                vPack_vec_final, vPack_vec_final_2, vPack_split, refOut, nItems,
                numRows, omega, h_rowDelimiters, filename, h_vec, Ntimes);
    }
}
#endif