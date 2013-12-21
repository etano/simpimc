#ifndef COMMUNICATION_H
#define COMMUNICATION_H

#include "../config.h"
#include <cstdint>


namespace COMM
{
#ifdef USE_MPI // Parallel version
  inline void Init (int argc, char **argv)
  {
    MPI_Init(&argc, &argv);
    int proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
  }

  inline void Finalize()
  {
    MPI_Finalize();
  }

  inline int WorldProc()
  {
    int proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
    return (proc);
  }

  inline void BarrierSync()
  {
    MPI_Barrier(MPI_COMM_WORLD);
  }

  // Template to retrieve traits of any MPI object
  template <class T>
  struct mpi_type_traits {
    static inline MPI_Datatype get_type(T&& val);
    static inline size_t get_size(T& val);
    static inline void* get_addr(T& val);
  };

  // Specialization of mpi_type_traits for primitive types
#define PRIMITIVE(Type, MpiType) \
        template<> \
        inline MPI_Datatype mpi_type_traits<Type>::get_type(Type&&) { return MpiType; } \
        template<> \
        inline size_t mpi_type_traits<Type>::get_size(Type&) { return 1; } \
        template<> \
        inline void* mpi_type_traits<Type>::get_addr(Type& val) { return &val; }
  PRIMITIVE(char, MPI::CHAR);
  PRIMITIVE(wchar_t, MPI::WCHAR);
  PRIMITIVE(short, MPI::SHORT);
  PRIMITIVE(int, MPI::INT);
  PRIMITIVE(long, MPI::LONG);
  PRIMITIVE(signed char, MPI::SIGNED_CHAR);
  PRIMITIVE(unsigned char, MPI::UNSIGNED_CHAR);
  PRIMITIVE(unsigned short, MPI::UNSIGNED_SHORT);
  PRIMITIVE(unsigned int, MPI::UNSIGNED);
  PRIMITIVE(unsigned long, MPI::UNSIGNED_LONG);
  PRIMITIVE(unsigned long long, MPI::UNSIGNED_LONG_LONG);
  PRIMITIVE(bool, MPI::BOOL);
#if PRECISION==double
  PRIMITIVE(RealType, MPI::DOUBLE);
  PRIMITIVE(ComplexType, MPI::DOUBLE_COMPLEX);
#elif PRECISION==single
  PRIMITIVE(RealType, MPI::FLOAT);
  PRIMITIVE(ComplexType, MPI::COMPLEX);
#endif

#undef PRIMITIVE

  // Specialization of mpi_type_traits for armadillo types
#define ARMATYPE(Type, ElemType, MpiType) \
        template<> \
        inline MPI_Datatype mpi_type_traits<Type>::get_type(Type&&) { return MpiType; } \
        template<> \
        inline size_t mpi_type_traits<Type>::get_size(Type& val) { return val.size(); } \
        template<> \
        inline void* mpi_type_traits<Type>::get_addr(Type& val) { return val.memptr(); }
  ARMATYPE(Imatrix, int, MPI::INT);
  ARMATYPE(Ivector, int, MPI::INT);
#if PRECISION==double
  ARMATYPE(Tmatrix, RealType, MPI::DOUBLE);
  ARMATYPE(Tvector, RealType, MPI::DOUBLE);
  ARMATYPE(Cmatrix, ComplexType, MPI::DOUBLE_COMPLEX);
  ARMATYPE(Cvector, ComplexType, MPI::DOUBLE_COMPLEX);
#elif PRECISION==single
  ARMATYPE(Tmatrix, RealType, MPI::FLOAT);
  ARMATYPE(Tvector, RealType, MPI::FLOAT);
  ARMATYPE(Cmatrix, ComplexType, MPI::COMPLEX);
  ARMATYPE(Cvector, ComplexType, MPI::COMPLEX);
#endif
#undef ARMATYPE

#else // Serial version
  inline void Init (int argc, char **argv) {}
  inline void Finalize () {}
  inline int WorldProc() {return (0);}
#endif
}


class CommStatusClass
{
public:
  int Source;
  int Tag;
  int Error;
  int Length;
};


class CommunicatorClass
{
public:

  CommunicatorClass()
  {
    SetWorld();
  }

#if USE_MPI // Parallel version
  MPI_Comm MPIComm;

  void SetWorld(); // Sets this communicator to be that of all the processes (i.e. MPI_WORLD)
  int MyProc();
  int NumProcs();
  string MyHost();
  void BarrierSync();
  void Split(int color, CommunicatorClass &newComm);
  void Subset(Imatrix &ranks, CommunicatorClass &newComm);

  /// Gets for MPI traits
  // Get type
  template <class T>
  inline MPI_Datatype GetMPIDatatype(T &val) { return COMM::mpi_type_traits<T>::get_type(val); }
  // Get address of first element
  template <class T>
  inline void* GetMPIAddr(T &val) { return COMM::mpi_type_traits<T>::get_addr(val); }
  // Get size of whole object
  template <class T>
  inline size_t GetMPISize(T &val) { return COMM::mpi_type_traits<T>::get_size(val); }

  // Send
  template<class T>
  inline int Send(int toProc, T &val)
  {
    return MPI_Send(GetMPIAddr(val), GetMPISize(val), GetMPIDatatype(val), toProc, 0, MPIComm);
  }

  // Receive
  template<class T>
  inline int Receive(int fromProc, T &val)
  {
    return MPI_Recv(GetMPIAddr(val), GetMPISize(val), GetMPIDatatype(val), fromProc, 0, MPIComm, MPI_STATUS_IGNORE);
  }

  // Sendrecv
  template<class T>
  inline int SendReceive (int fromProc, T &fromBuff, int toProc, T &toBuff)
  {
    return MPI_Sendrecv(GetMPIAddr(fromBuff), GetMPISize(fromBuff), GetMPIDatatype(fromBuff), fromProc, 1,
                        GetMPIAddr(toBuff), GetMPISize(toBuff), GetMPIDatatype(toBuff), toProc, 1, MPIComm, MPI_STATUS_IGNORE);
  }

  // Broadcast
  template<class T>
  inline int Broadcast(int fromProc, T &val)
  {
    return MPI_Bcast(GetMPIAddr(val), GetMPISize(val), GetMPIDatatype(val), fromProc, MPIComm);
  }

  // Reduce
  template<class T>
  inline int Reduce(int toProc, T &fromBuff, T &toBuff, MPI_Op Op)
  {
    return MPI_Reduce(GetMPIAddr(fromBuff), GetMPIAddr(toBuff), GetMPISize(fromBuff), GetMPIDatatype(fromBuff),
                      Op, toProc, MPIComm);
  }

  // AllReduce
  template<class T>
  inline int AllReduce(T &fromBuff, T &toBuff, MPI_Op Op)
  {
    return MPI_Allreduce(GetMPIAddr(fromBuff), GetMPIAddr(toBuff), GetMPISize(fromBuff), GetMPIDatatype(fromBuff),
                      Op, MPIComm);
  }

  // Sum
  template<class T>
  inline int Sum(int toProc, T &fromBuff, T &toBuff)
  {
    return Reduce(toProc,fromBuff,toBuff,MPI_SUM);
  }

  // AllSum
  template<class T>
  inline int AllSum(T &fromBuff, T &toBuff)
  {
    return AllReduce(fromBuff,toBuff,MPI_SUM);
  }

  // Product
  template<class T>
  inline int Product(int toProc, T &fromBuff, T &toBuff)
  {
    return Reduce(toProc,fromBuff,toBuff,MPI_PROD);
  }

  // AllProduct
  template<class T>
  inline int AllProduct(T &fromBuff, T &toBuff)
  {
    return AllReduce(fromBuff,toBuff,MPI_PROD);
  }

  // Gather
  template<class T>
  inline int Gather(int toProc, T &fromBuff, T &toBuff)
  {
    return MPI_Gather(GetMPIAddr(fromBuff), GetMPISize(fromBuff), GetMPIDatatype(fromBuff),
                      GetMPIAddr(toBuff), GetMPISize(fromBuff), GetMPIDatatype(toBuff), toProc, MPIComm);
  }

  // Gatherv
  template<class T>
  inline int Gatherv(int toProc, T &fromBuff, T &toBuff, int* recvCounts, int* displacements)
  {
     return MPI_Gatherv(GetMPIAddr(fromBuff), GetMPISize(fromBuff), GetMPIDatatype(fromBuff),
                        GetMPIAddr(toBuff), recvCounts, displacements, GetMPIDatatype(toBuff), toProc, MPIComm);
  }

  // GatherCols
  template<class T>
  int GatherCols(int toProc, T &fromBuff, T &toBuff)
  {
    int nProcs = NumProcs();
    int myProc = MyProc();
    int rows = fromBuff.n_rows;
    int cols = fromBuff.n_cols;
    int displacements[nProcs];
    int recvCounts[nProcs];
    int currCol = 0;
    for (int proc=0; proc<nProcs; proc++) {
      int procCols = cols/nProcs + ((cols%nProcs)>proc);
      displacements[proc] = rows*currCol;
      recvCounts[proc] = rows*procCols;
      if (proc == myProc) {
        toBuff.set_size(rows,procCols);
      }
      currCol += procCols;
    }
    return Gather(toProc, fromBuff, toBuff, recvCounts, displacements);
  }

  // AllGather
  template<class T>
  inline int AllGather(T &fromBuff, T &toBuff)
  {
    return MPI_Gather(GetMPIAddr(fromBuff), GetMPISize(fromBuff), GetMPIDatatype(fromBuff),
                      GetMPIAddr(toBuff), GetMPISize(fromBuff), GetMPIDatatype(toBuff), MPIComm);
  }

  // Scatter
  template<class T>
  inline int Scatter(int fromProc, T &fromBuff, T &toBuff)
  {
    return MPI_Scatter(GetMPIAddr(fromBuff), GetMPISize(fromBuff), GetMPIDatatype(fromBuff),
                       GetMPIAddr(toBuff), GetMPISize(toBuff), GetMPIDatatype(toBuff), fromProc, MPIComm);
  }

  // Scatter
  template<class T>
  inline int Scatterv(int fromProc, T &fromBuff, T &toBuff, int* sendCounts, int* displacements)
  {
     return MPI_Scatterv(GetMPIAddr(fromBuff), sendCounts, displacements, GetMPIDatatype(fromBuff),
                         GetMPIAddr(toBuff), GetMPISize(toBuff), GetMPIDatatype(toBuff), fromProc, MPIComm);
  }

  // ScatterCols
  template<class T>
  int ScatterCols(int fromProc, T &fromBuff, T &toBuff)
  {
    int nProcs = NumProcs();
    int myProc = MyProc();
    int rows = fromBuff.n_rows;
    int cols = fromBuff.n_cols;
    int displacements[nProcs];
    int sendCounts[nProcs];
    int currCol = 0;
    for (int proc=0; proc<nProcs; proc++) {
      int procCols = cols/nProcs + ((cols%nProcs)>proc);
      displacements[proc] = rows*currCol;
      sendCounts[proc] = rows*procCols;
      if (proc == myProc)
        toBuff.set_size(rows,procCols);
      currCol += procCols;
    }
    return Scatterv(fromProc, fromBuff, toBuff, sendCounts, displacements);
  }

  // AllGatherCols
  int AllGatherCols (Tmatrix &buff);

#else   // Serial version
  inline void SetWorld(){}
  inline int MyProc() {return 0;}
  inline int NumProcs() {return 1;}
  inline string MyHost() {return "0";}
  inline void BarrierSync() {}
  inline void Split(int color, CommunicatorClass &newComm) {}
  inline void Subset(Imatrix ranks, CommunicatorClass &newComm)
  {
    if (ranks.size() != 1) {
      cerr << "Serial verion of code does not support nontrivial subsets. Exiting.\n";
      abort();
    }
  }

  template<class T>
  inline int Send(int toProc, T &val) {}
  template<class T>
  inline int Receive(int fromProc, T &val) {}
  template<class T>
  inline int SendReceive (int fromProc, T &fromBuff, int toProc, T &toBuff) {toBuff = fromBuff;}
  template<class T>
  inline int Broadcast(int fromProc, T &val) {}
  template<class T>
  inline int Sum(int toProc, T &fromBuff, T &toBuff) {toBuff = fromBuff;}
  template<class T>
  inline int AllSum(T &fromBuff, T &toBuff) {toBuff = fromBuff;}
  template<class T>
  inline int Product(int toProc, T &fromBuff, T &toBuff) {toBuff = fromBuff;}
  template<class T>
  inline int AllProduct(T &fromBuff, T &toBuff) {toBuff = fromBuff;}
  template<class T>
  inline int Gather(int toProc, T &fromBuff, T &toBuff) {toBuff = fromBuff;}
  template<class T>
  inline int Gatherv(int toProc, T &fromBuff, T &toBuff, int* recvCounts, int* displacements) {toBuff = fromBuff;}
  template<class T>
  inline int GatherCols(int toProc, T &fromBuff, T &toBuff) {toBuff = fromBuff;}
  template<class T>
  inline int AllGather(T &fromBuff, T &toBuff) {toBuff = fromBuff;}
  template<class T>
  inline int Scatter(int fromProc, T &fromBuff, T &toBuff) {toBuff = fromBuff;}
  template<class T>
  inline int Scatterv(int fromProc, T &fromBuff, T &toBuff, int* sendCounts, int* displacements) {toBuff = fromBuff;}
  template<class T>
  inline int ScatterCols(int fromProc, T &fromBuff, T &toBuff) {toBuff = fromBuff;}
  inline int AllGatherCols (Tmatrix &buff) {}
#endif
};

#endif

