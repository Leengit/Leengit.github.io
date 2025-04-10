// Optimal Checkpointing for Backtrace.
// 
// Author: Lee A. Newberg.
// Wadsworth Center, New York State Department of Health, Albany NY USA
// Dept. of Computer Science, Rensselaer Polytechnic Institute, Troy NY USA

#include <iostream>

// Warning: the BIGINT type, must be able to handle values a little
// more than N*M*M, otherwise this software will produce unpredictable
// results.  For instance, N=10^6, M=10^3 is too large for a 4-byte
// integer.
typedef long long BIGINT;

//
// Prototypes
//

// main: a sample program that performs a backtrace.
int main(int argc, char *argv[]);

// backtrace: the workhorse routine that performs the backtrace.  The
// initial value for each of Mckpt and Nckpt should be one less than
// the lowest array index, e.g., the initial value should be -1 for
// typical C and C++ implementations.  The initial value for L and
// Nopt should be as computed from the findLevel routine.
void backtrace
(  BIGINT Mckpt,                // previous checkpoint memory location
   BIGINT M,                    // M unused memory locations
   BIGINT Nckpt,                // previous checkpoint stage
   BIGINT N,                    // backtrace through N stages
   BIGINT L,                    // applicable special level
   BIGINT Nopt,                 // applicable special stage
   // advance, available, and p are for user callbacks:
   void (*advance)(BIGINT Mfrom, BIGINT Mto, BIGINT Nto, void *p),
   void (*available)(BIGINT M, BIGINT N, void *p),
   void *p);                    // User-supplied global information

// findLevel: compute level L and Nopt(M,L) from initial values for M
// and N.
void findLevel(BIGINT M, BIGINT N, BIGINT &L, BIGINT &Nopt);

// advance: a user callback function for computing a stage from its
// immediate predecessor.  The supplied version of this function
// merely prints the task it was asked to perform.  Normally, the
// argument p would point to a structure that included a pointer to
// the memory available for stage computations, and the task would be
// peformed in that memory.  When Mfrom is one lower than the lowest
// index for a memory location (e.g., when it is -1, which is one
// lower than 0, the first array index in C and C++ implementations)
// then this callback function should compute Nto from
// initial/boundary conditions.
void myAdvance(BIGINT Mfrom, BIGINT Mto, BIGINT Nto, void *p);

// available: a user callback for when a stage becomes available
// during the backtrace.  The supplied version of this function merely
// prints the task it was asked to perform.  Normally, the argument p
// would point to a structure that included a pointer to the memory
// available for stage computations, and the available stage would be
// read from that memory.
void myAvailable(BIGINT Mavailable, BIGINT Navailable, void *p);

// computeTopt: a method that computes the required number of calls to
// the advance callback routine.  This routine is not used by the
// implementation of backtrace, but is for the user's information
// only.
BIGINT computeTopt(BIGINT M, BIGINT N, BIGINT L, BIGINT Nopt);

//
// Function definitions
//

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        std::cout << std::flush;
        std::cerr << argv[0]
                  << ": Please supply two command line arguments, M and N."
                  << std::endl << std::flush;
        exit(-1);
    }
    BIGINT Mckpt = -1;
    BIGINT M = atol(argv[1]);
    BIGINT Nckpt = -1;
    BIGINT N = atol(argv[2]);
    BIGINT Nopt, L;
    findLevel(M, N, L, Nopt);
    int p[1];
    p[0] = 0;
    std::cout << "Topt: "
              << computeTopt(M, N, L, Nopt)
              << " calls to the advance callback routine will be required."
              << std::endl;
    std::cout << "Calling backtrace("
              << "Mckpt = " << Mckpt << ", "
              << "M = " << M << ", "
              << "Nckpt = " << Nckpt << ", "
              << "N = " << N << ", "
              << "L = " << L << ", "
              << "Nopt = " << Nopt << ", "
              << std::endl << "\t"
              << "advance = &myAdvance, "
              << "available = &myAvailable, "
              << "p = " << p << ")."
              << std::endl;
    backtrace(Mckpt, M, Nckpt, N, L, Nopt, &myAdvance, &myAvailable, p);
    std::cout << "Verify Topt: "
              << p[0]
              << " calls to the advance callback routine were required."
              << std::endl;
}

void findLevel(BIGINT M, BIGINT N, BIGINT &L, BIGINT &Nopt)
{
    // Warning: this routine should start with a better guess for L if
    // it is likely to be used for M <= 3 and very large N.
    if (!(M >= 2 && N >= 1))
    {
        std::cout << std::flush;
        std::cerr << "Bad inputs (M, N) = (" << M << ", " << N << ")."
                  << std::endl
                  << "M should be at least 2.  N should be at least 1."
                  << std::endl << std::flush;
        exit(-2);
    }
    L = 0;
    Nopt = 0;
    BIGINT Lnew = 1;
    BIGINT Noptnew = M;
    while (Noptnew <= N)
    {
        L = Lnew;
        Nopt = Noptnew;
        Lnew = L + 1;
        Noptnew = ((Nopt+1)*(M+L-1)*(M+L*2+1)) / ((M+L*2-1)*(L+1)) - 1;
    }
}

void myAdvance(BIGINT Mfrom, BIGINT Mto, BIGINT Nto, void *p)
{
    BIGINT Nfrom = Nto - 1;
    std::cout << "Advance: From stage " << Nfrom
              << " (location " << Mfrom
              << ") computing stage " << Nto
              << " (location " << Mto
              << ")."
              << std::endl;
    // Use the user-supplied pointer.
    ++*reinterpret_cast<int *>(p); // Count calls to myAdvance
}

void myAvailable(BIGINT Mavailable, BIGINT Navailable, void *p)
{
    std::cout << "Available: Backtracing to stage " << Navailable
              << " (location " << Mavailable
              << ")."
              << std::endl;
}

BIGINT computeTopt(BIGINT M, BIGINT N, BIGINT L, BIGINT Nopt)
{
    return ((Nopt+1) * (M*M+M*L*2-M*2-L*2+2) * L) / ((M+L*2-1) * M)
        + (L+1) * (N-Nopt);
}

void
backtrace(BIGINT Mckpt,
          BIGINT M,
          BIGINT Nckpt,
          BIGINT N,
          BIGINT L,
          BIGINT Nopt,
          void (*advance)(BIGINT Mfrom, BIGINT Mto, BIGINT Nto, void *p),
          void (*available)(BIGINT M, BIGINT N, void *p),
          void *p)
{
    while (N > M)
    {
        BIGINT Noptm;           // Nopt(M-1,L)
        BIGINT Noptl;           // Nopt(M,L-1)
        BIGINT C;               // Nckpt + C is the next stage to checkpoint
        Noptm = ((Nopt+1)*(M-1)*(M+L*2-2)) / ((M+L*2-1)*(M+L-2)) - 1;
        Noptl = ((Nopt+1)*L*(M+L*2-3)) / ((M+L*2-1)*(M+L-2)) - 1;
        C = (Nopt+1 < N-Noptm ? Nopt+1 : N-Noptm);
        // Compute stage Nckpt + C from stage Nckpt, by alternate use of
        // two memory locations Mckpt+1 and Mckpt +2
        (*advance)(Mckpt, Mckpt+1+((C-1)%2), Nckpt+1, p);
        for (BIGINT i = 2; i <= C; ++i)
            (*advance)(Mckpt+2-((C-i)%2), Mckpt+1+((C-i)%2),
                       Nckpt+i, p);
        // Backtrace through stages Nckpt+C+1, ..., Nckpt+N
        backtrace(Mckpt+1,M-1,Nckpt+C,N-C,L,Noptm,advance,available,p);
        // Present stage Nckpt + C to the user
        (*available)(Mckpt+1, Nckpt + C, p);
        // Backtrace through stages Nckpt+1, ... Nckpt+C-1 via tail recursion
        N = C - 1;
        L = L - 1;
        Nopt = Noptl;
    }
    // Handle ample memory case, N <= M
    for (BIGINT i = 1; i <= N; ++i)
        (*advance)(Mckpt+i-1, Mckpt+i, Nckpt+i, p);
    for (BIGINT i = N; i >= 1; --i)
        (*available)(Mckpt+i, Nckpt+i, p);
}
