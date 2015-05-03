#ifndef LAPACK
#define LAPACK

//#define F77NAME(x) x##_
//void F77NAME(dgetrf)(const long int* m, const long int* n, double* a, const long int* lda, long int* ipiv, long int* info);
//void F77NAME(dgetrs)(const char* trans, const long int* n, const long int* nrhs, const double* a, const long int* lda, const long int* ipiv, double* b, const long int* ldb, long int* info);

extern void dgetrf_(const long int* m, const long int* n, double* a, const long int* lda, long int* ipiv, long int* info);
extern void dgetrs_(const char* trans, const long int* n, const long int* nrhs, const double* a, const long int* lda, const long int* ipiv, double* b, const long int* ldb, long int* info);
extern void dgesv_(const long int* n, const long int* nrhs, double* a, const long int* lda, long int* ipiv, double* b, const long int* ldb, long int* info);

#endif //LAPACK
