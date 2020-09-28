#ifndef _LBLAS_
#define _LBLAS_
void dsygvd_(const int* ITYPE, const char* JOBZ, const char* UPLO, const int* N,
    double* A, const int* LDA, double* B, const int* LDB, double* W, double* WORK, const int* LWORK,
    int* IWORK, const int* LIWORK, int* INFO)
{

}
int Eigensystem(int N, function<double>& Energy, const function2D<double>& Olap, function2D<double>& F)
{
    static function2D<double> tOlap(N, N);
    int lwork = 1 + 6 * N + 2 * pow(N, 2) + 10;
    static function1D<double> work(lwork);
    int liwork = 3 + 5 * N + 10;
    static function1D<int> iwork(liwork);
    int itype = 1;
    int info = 0;
    tOlap = Olap;
    int lF = F.SizeNd();
    int lO = Olap.SizeNd();
    dsygvd_(&itype, "V", "U", &N, F.MemPt(), &lF, tOlap.MemPt(), &lO, Energy.MemPt(), work.MemPt(), &lwork, iwork.MemPt(), &liwork, &info);
    if (info)
    {
        cerr << "Solving eigenvalue problem is unsuccessful " << info << endl;
    }
    return info;
}
#endif //_LBLAS_