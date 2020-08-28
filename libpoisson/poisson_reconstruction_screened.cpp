#include "poisson_reconstruction_screened.h"

#pragma warning(push)
#pragma warning(disable:4100)
#pragma warning(disable:4127)
#pragma warning(disable:4189)
#pragma warning(disable:4244)
#pragma warning(disable:4701)
#pragma warning(disable:4702)
#pragma warning(disable:4706)
//#define _OPENMP

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#ifdef _WIN32
#include <Windows.h>
#include <Psapi.h>
#endif

#include "poisson_screened/MyTime.h"
#include "poisson_screened/MarchingCubes.h"
#include "poisson_screened/Octree.h"
#include "poisson_screened/SparseMatrix.h"
#include "poisson_screened/CmdLineParser.h"
#include "poisson_screened/PPolynomial.h"
#include "poisson_screened/PlyVertexMini.h"
#include "poisson_screened/MemoryUsage.h"
#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP
#include "poisson_screened/MultiGridOctreeData.h"

#include <iostream>
#include <cassert>

namespace
  {
  ///////////////////////////////////////////////////////////////////////////////

  Point3D<float> Convert(const jtk::vec3<float>& pt)
    {
    return Point3D<float>(pt[0], pt[1], pt[2]);
    }

  ///////////////////////////////////////////////////////////////////////////////

  class MyPointStream : public OrientedPointStream < float >
    {
    public:
      MyPointStream(const std::vector<jtk::vec3<float>>& pts, const std::vector<jtk::vec3<float>>& normals) : p_pts(&pts), p_normals(&normals) {}
      void reset()
        {
        _current = 0;
        }
      virtual bool nextPoint(OrientedPoint3D< float >& p)
        {
        if (_current >= p_pts->size())
          return false;
        p.p = Convert((*p_pts)[_current]);
        p.n = Convert((*p_normals)[_current]);
        ++_current;
        return true;
        }
    private:
      const std::vector<jtk::vec3<float>>* p_pts;
      const std::vector<jtk::vec3<float>>* p_normals;
      int _current;
    };

  ///////////////////////////////////////////////////////////////////////////////

  class MyColoredPointStream : public OrientedPointStreamWithData < float, Point3D<float> >
    {
    public:
      MyColoredPointStream(const std::vector<jtk::vec3<float>>& pts, const std::vector<jtk::vec3<float>>& normals, const std::vector<uint32_t>& colors)
        : p_pts(&pts), p_normals(&normals), p_colors(&colors) {}
      void reset()
        {
        _current = 0;
        }
      virtual bool nextPoint(OrientedPoint3D< float >& p, Point3D<float>& d)
        {
        if (_current >= p_pts->size())
          return false;
        p.p = Convert((*p_pts)[_current]);
        p.n = Convert((*p_normals)[_current]);
        uint32_t clr = (*p_colors)[_current];
        uint32_t red = clr & 255;
        uint32_t green = (clr >> 8) & 255;
        uint32_t blue = (clr >> 16) & 255;
        d[0] = red;
        d[1] = green;
        d[2] = blue;
        ++_current;
        return true;
        }
    private:
      const std::vector<jtk::vec3<float>>* p_pts;
      const std::vector<jtk::vec3<float>>* p_normals;
      const std::vector<uint32_t>* p_colors;
      int _current;
    };

  ///////////////////////////////////////////////////////////////////////////////

  void DumpOutput(const char* format, ...)
    {
    char buf[4096];
    va_list marker;
    va_start(marker, format);

    vsprintf(buf, format, marker);
    va_end(marker);

    std::cout << buf;
    }
  void DumpOutput2(std::vector< char* >& comments, const char* format, ...)
    {
    char buf[4096];
    va_list marker;
    va_start(marker, format);

    vsprintf(buf, format, marker);
    va_end(marker);
    std::cout << buf;
    }


  template< class Real>
  XForm4x4<Real> GetPointStreamScale(const std::vector<jtk::vec3<Real>>& vertices, Real expFact)
    {
    auto bb = jtk::bounding_volume_3d<float>(vertices.begin(), vertices.end());

    Real scale = std::max<Real>(bb.max[0] - bb.min[0], std::max<Real>(bb.max[1] - bb.min[1], bb.max[2] - bb.min[2]));
    scale *= expFact;
    //Real scale = bb.Dim()[bb.MaxDim()] * expFact;
    auto center = jtk::center(bb);
    for (int i = 0; i < 3; i++) center[i] -= scale / 2;
    XForm4x4< Real > tXForm = XForm4x4< Real >::Identity(), sXForm = XForm4x4< Real >::Identity();
    for (int i = 0; i < 3; i++) sXForm(i, i) = (Real)(1. / scale), tXForm(3, i) = -center[i];
    return sXForm * tXForm;
    }
  }




void poisson_reconstruction_screened(std::vector<jtk::vec3<float>>& vertices, std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<float>>& pts, const std::vector<jtk::vec3<float>>& normals, const poisson_reconstruction_screened_parameters& par)
  {
  int MaxDepthVal = 8;
  int MaxSolveDepthVal = -1;
  int KernelDepthVal = -1;
  int MinDepthVal = 0;
  int FullDepthVal = 5;
  float SamplesPerNodeVal = 1.5f;
  float ScaleVal = 1.1f;
  bool ConfidenceFlag = false;
  bool CleanFlag = false;
  bool DensityFlag = false;
  float PointWeightVal = 4.f;
  int AdaptiveExponentVal = 1;
  int BoundaryTypeVal = 1;
  bool CompleteFlag = false;
  bool NonManifoldFlag = false;
  bool ShowResidualFlag = false;
  int CGDepthVal = 0;
  int ItersVal = 8;
  float CSSolverAccuracyVal = 1e-3f;

  bool VerboseFlag = true;
  int ThreadsVal = omp_get_num_procs();
  bool LinearFitFlag = false;
  float LowResIterMultiplierVal = 1.f;
  float ColorVal = 16.f;

#define Degree 2

  typedef float Real;

  typedef typename Octree< Real >::template DensityEstimator< WEIGHT_DEGREE > DensityEstimator;
  typedef typename Octree< Real >::template InterpolationInfo< false > InterpolationInfo;
  //typedef OrientedPointStreamWithData< Real, Point3D< Real > > PointStreamWithData;
  typedef TransformedOrientedPointStream< Real > XPointStream;
  Reset< Real >();
  std::vector< char* > comments;

  XForm4x4< Real > xForm = GetPointStreamScale<Real>(pts, ScaleVal);
  XForm4x4< Real > iXForm = xForm.inverse();
  DumpOutput2(comments, "Running Screened Poisson Reconstruction (Version 9.0)\n");
  double startTime = Time();

  OctNode< TreeNodeData >::SetAllocator(MEMORY_ALLOCATOR_BLOCK_SIZE);
  Octree< Real > tree;
  //OctreeProfiler< Real > profiler(tree);
  tree.threads = ThreadsVal;
  if (MaxSolveDepthVal < 0) MaxSolveDepthVal = MaxDepthVal;




  //	int kernelDepth = KernelDepth.set ? KernelDepth.value : Depth.value-2;
  if (KernelDepthVal < 0) KernelDepthVal = MaxDepthVal - 2;
  if (KernelDepthVal > MaxDepthVal)
    {
    printf("kernelDepth cannot be greateer Depth.value\n");
    return;
    }

  int pointCount;

  Real pointWeightSum;
  std::vector< typename Octree< Real >::PointSample >* samples = new std::vector< typename Octree< Real >::PointSample >();
  std::vector< ProjectiveData< Point3D< Real >, Real > >* sampleData = NULL;
  DensityEstimator* density = NULL;
  SparseNodeData< Point3D< Real >, NORMAL_DEGREE >* normalInfo = NULL;
  Real targetValue = (Real)0.5;

  // Read in the samples (and color data)
  {

  sampleData = nullptr;// new std::vector< ProjectiveData< Point3D< Real >, Real > >();
  MyPointStream my_pointStream(pts, normals);
  my_pointStream.reset();
  XPointStream _pointStream(xForm, my_pointStream);
  pointCount = tree.template init< Point3D< Real > >(_pointStream, MaxDepthVal, ConfidenceFlag, *samples, sampleData);

#pragma omp parallel for num_threads( ThreadsVal )
  for (int i = 0; i < (int)samples->size(); i++) (*samples)[i].sample.data.n *= (Real)-1;

  DumpOutput("Input Points / Samples: %d / %d\n", pointCount, samples->size());
  //profiler.dumpOutput2(comments, "# Read input into tree:");
  }

  DenseNodeData< Real, Degree > solution;

  {
  DenseNodeData< Real, Degree > constraints;
  InterpolationInfo* iInfo = NULL;
  int solveDepth = MaxSolveDepthVal;

  tree.resetNodeIndices();

  // Get the kernel density estimator [If discarding, compute anew. Otherwise, compute once.]
  {
  //profiler.start();
  density = tree.template setDensityEstimator< WEIGHT_DEGREE >(*samples, KernelDepthVal, SamplesPerNodeVal);
  //profiler.dumpOutput2(comments, "#   Got kernel density:");
  }

  // Transform the Hermite samples into a vector field [If discarding, compute anew. Otherwise, compute once.]
  {
  //profiler.start();
  normalInfo = new SparseNodeData< Point3D< Real >, NORMAL_DEGREE >();
  *normalInfo = tree.template setNormalField< NORMAL_DEGREE >(*samples, *density, pointWeightSum, BOUNDARY_NEUMANN == BOUNDARY_NEUMANN);
  //profiler.dumpOutput2(comments, "#     Got normal field:");
  }

  if (!DensityFlag) delete density, density = NULL;

  // Trim the tree and prepare for multigrid
  {
  //profiler.start();
  std::vector< int > indexMap;

  constexpr int MAX_DEGREE = NORMAL_DEGREE > Degree ? NORMAL_DEGREE : Degree;
  tree.template inalizeForBroodedMultigrid< MAX_DEGREE, Degree, BOUNDARY_NEUMANN >(FullDepthVal, typename Octree< Real >::template HasNormalDataFunctor< NORMAL_DEGREE >(*normalInfo), &indexMap);

  if (normalInfo) normalInfo->remapIndices(indexMap);
  if (density) density->remapIndices(indexMap);
  //profiler.dumpOutput2(comments, "#       Finalized tree:");
  }

  // Add the FEM constraints
  {
  //profiler.start();
  constraints = tree.template initDenseNodeData< Degree >();
  tree.template addFEMConstraints< Degree, BOUNDARY_NEUMANN, NORMAL_DEGREE, BOUNDARY_NEUMANN >(FEMVFConstraintFunctor< NORMAL_DEGREE, BOUNDARY_NEUMANN, Degree, BOUNDARY_NEUMANN >(1., 0.), *normalInfo, constraints, solveDepth);
  //profiler.dumpOutput2(comments, "#  Set FEM constraints:");
  }

  // Free up the normal info [If we don't need it for subseequent iterations.]
  delete normalInfo, normalInfo = NULL;

  // Add the interpolation constraints
  if (PointWeightVal > 0)
    {
    //profiler.start();
    iInfo = new InterpolationInfo(tree, *samples, targetValue, AdaptiveExponentVal, (Real)PointWeightVal * pointWeightSum, (Real)0);
    tree.template addInterpolationConstraints< Degree, BOUNDARY_NEUMANN >(*iInfo, constraints, solveDepth);
    //profiler.dumpOutput2(comments, "#Set point constraints:");
    }

  DumpOutput("Leaf Nodes / Active Nodes / Ghost Nodes: %d / %d / %d\n", (int)tree.leaves(), (int)tree.nodes(), (int)tree.ghostNodes());
  DumpOutput("Memory Usage: %.3f MB\n", float(MemoryInfo::Usage()) / (1 << 20));

  // Solve the linear system
  {
  //profiler.start();
  typename Octree< Real >::SolverInfo solverInfo;
  solverInfo.cgDepth = CGDepthVal, solverInfo.iters = ItersVal, solverInfo.cgAccuracy = CSSolverAccuracyVal, solverInfo.verbose = VerboseFlag, solverInfo.showResidual = ShowResidualFlag, solverInfo.lowResIterMultiplier = std::max< double >(1., LowResIterMultiplierVal);
  solution = tree.template solveSystem< Degree, BOUNDARY_NEUMANN >(FEMSystemFunctor< Degree, BOUNDARY_NEUMANN >(0, 1., 0), iInfo, constraints, solveDepth, solverInfo);
  //profiler.dumpOutput2(comments, "# Linear system solved:");
  if (iInfo) delete iInfo, iInfo = NULL;
  }
  }

  typedef PlyVertex< float > Vertex;

  CoredFileMeshData< Vertex > mesh;

  {
  //profiler.start();
  double valueSum = 0, weightSum = 0;
  typename Octree< Real >::template MultiThreadedEvaluator< Degree, BOUNDARY_NEUMANN > evaluator(&tree, solution, ThreadsVal);
#pragma omp parallel for num_threads( ThreadsVal ) reduction( + : valueSum , weightSum )
  for (int j = 0; j < samples->size(); j++)
    {
    ProjectiveData< OrientedPoint3D< Real >, Real >& sample = (*samples)[j].sample;
    Real w = sample.weight;
    if (w > 0) weightSum += w, valueSum += evaluator.value(sample.data.p / sample.weight, omp_get_thread_num(), (*samples)[j].node) * w;
    }
  Real isoValue = (Real)(valueSum / weightSum);
  //		if( samples ) delete samples , samples = NULL;
  //profiler.dumpOutput("Got average:");
  DumpOutput("Iso-Value: %e\n", isoValue);

  //profiler.start();
  SparseNodeData< ProjectiveData< Point3D< Real >, Real >, DATA_DEGREE >* colorData = NULL;
  if (sampleData)
    {
    colorData = new SparseNodeData< ProjectiveData< Point3D< Real >, Real >, DATA_DEGREE >();
    *colorData = tree.template setDataField< DATA_DEGREE, false >(*samples, *sampleData, (DensityEstimator*)NULL);
    delete sampleData, sampleData = NULL;
    for (const OctNode< TreeNodeData >* n = tree.tree().nextNode(); n; n = tree.tree().nextNode(n))
      {
      ProjectiveData< Point3D< Real >, Real >* clr = (*colorData)(n);
      if (clr)
        (*clr) *= (Real)pow(ColorVal, tree.depth(n));
      }
    }
  tree.template getMCIsoSurface< Degree, BOUNDARY_NEUMANN, WEIGHT_DEGREE, DATA_DEGREE >(density, colorData, solution, isoValue, mesh, !LinearFitFlag, !NonManifoldFlag, false /*PolygonMesh.set*/);
  DumpOutput("Vertices / Polygons: %d / %d\n", mesh.outOfCorePointCount() + mesh.inCorePoints.size(), mesh.polygonCount());
  //profiler.dumpOutput2(comments, "#        Got triangles:");
  }

  //        FreePointer( solution );


  mesh.resetIterator();
  int nr_vertices = int(mesh.outOfCorePointCount() + mesh.inCorePoints.size());
  int nr_faces = mesh.polygonCount();

  vertices.clear();
  triangles.clear();

  for (int i = 0; i<int(mesh.inCorePoints.size()); ++i)
    {
    PlyVertex< float > vertex = mesh.inCorePoints[i];

    Point3D<Real> pp = iXForm * vertex.point;

    vertices.push_back(jtk::vec3<float>(pp[0], pp[1], pp[2]));
    }
  for (int i = 0; i < mesh.outOfCorePointCount(); ++i)
    {
    PlyVertex< float > vertex;
    mesh.nextOutOfCorePoint(vertex);
    Point3D<Real> pp = iXForm * vertex.point;

    vertices.push_back(jtk::vec3<float>(pp[0], pp[1], pp[2]));
    }
  std::vector< CoredVertexIndex > polygon;
  while (mesh.nextPolygon(polygon))
    {
    assert(polygon.size() == 3);
    int indV[3];
    for (int i = 0; i<int(polygon.size()); i++)
      {
      if (polygon[i].inCore) indV[i] = polygon[i].idx;
      else                    indV[i] = polygon[i].idx + int(mesh.inCorePoints.size());
      }
    jtk::vec3<uint32_t> tria(indV[0], indV[1], indV[2]);
    triangles.push_back(tria);
    }

  if (density) delete density, density = NULL;
  DumpOutput2(comments, "#          Total Solve: %9.1f (s), %9.1f (MB)\n", Time() - startTime, tree.maxMemoryUsage());


  }



void poisson_reconstruction_screened(std::vector<jtk::vec3<float>>& vertices, std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<jtk::vec3<float>>& vertex_colors, const std::vector<jtk::vec3<float>>& pts, const std::vector<jtk::vec3<float>>& normals, const std::vector<uint32_t>& colors, const poisson_reconstruction_screened_parameters& par)
  {
  int MaxDepthVal = 8;
  int MaxSolveDepthVal = -1;
  int KernelDepthVal = -1;
  int MinDepthVal = 0;
  int FullDepthVal = 5;
  float SamplesPerNodeVal = 1.5f;
  float ScaleVal = 1.1f;
  bool ConfidenceFlag = false;
  bool CleanFlag = false;
  bool DensityFlag = false;
  float PointWeightVal = 4.f;
  int AdaptiveExponentVal = 1;
  int BoundaryTypeVal = 1;
  bool CompleteFlag = false;
  bool NonManifoldFlag = false;
  bool ShowResidualFlag = false;
  int CGDepthVal = 0;
  int ItersVal = 8;
  float CSSolverAccuracyVal = 1e-3f;

  bool VerboseFlag = true;
  int ThreadsVal = omp_get_num_procs();
  bool LinearFitFlag = false;
  float LowResIterMultiplierVal = 1.f;
  float ColorVal = 16.f;

#define Degree 2

  typedef float Real;

  typedef typename Octree< Real >::template DensityEstimator< WEIGHT_DEGREE > DensityEstimator;
  typedef typename Octree< Real >::template InterpolationInfo< false > InterpolationInfo;
  //typedef OrientedPointStreamWithData< Real, Point3D< Real > > PointStreamWithData;
  typedef TransformedOrientedPointStreamWithData< Real, Point3D<Real> > XPointStreamWithData;
  Reset< Real >();
  std::vector< char* > comments;

  XForm4x4< Real > xForm = GetPointStreamScale<Real>(pts, ScaleVal);
  XForm4x4< Real > iXForm = xForm.inverse();
  DumpOutput2(comments, "Running Screened Poisson Reconstruction (Version 9.0)\n");
  double startTime = Time();

  OctNode< TreeNodeData >::SetAllocator(MEMORY_ALLOCATOR_BLOCK_SIZE);
  Octree< Real > tree;
  //OctreeProfiler< Real > profiler(tree);
  tree.threads = ThreadsVal;
  if (MaxSolveDepthVal < 0) MaxSolveDepthVal = MaxDepthVal;




  //	int kernelDepth = KernelDepth.set ? KernelDepth.value : Depth.value-2;
  if (KernelDepthVal < 0) KernelDepthVal = MaxDepthVal - 2;
  if (KernelDepthVal > MaxDepthVal)
    {
    printf("kernelDepth cannot be greateer Depth.value\n");
    return;
    }

  int pointCount;

  Real pointWeightSum;
  std::vector< typename Octree< Real >::PointSample >* samples = new std::vector< typename Octree< Real >::PointSample >();
  std::vector< ProjectiveData< Point3D< Real >, Real > >* sampleData = NULL;
  DensityEstimator* density = NULL;
  SparseNodeData< Point3D< Real >, NORMAL_DEGREE >* normalInfo = NULL;
  Real targetValue = (Real)0.5;

  // Read in the samples (and color data)
  {

  sampleData = new std::vector< ProjectiveData< Point3D< Real >, Real > >();
  MyColoredPointStream my_pointStream(pts, normals, colors);
  my_pointStream.reset();
  XPointStreamWithData _pointStream(xForm, my_pointStream);
  pointCount = tree.template init< Point3D< Real > >(_pointStream, MaxDepthVal, ConfidenceFlag, *samples, sampleData);

#pragma omp parallel for num_threads( ThreadsVal )
  for (int i = 0; i < (int)samples->size(); i++) (*samples)[i].sample.data.n *= (Real)-1;

  DumpOutput("Input Points / Samples: %d / %d\n", pointCount, samples->size());
  //profiler.dumpOutput2(comments, "# Read input into tree:");
  }

  DenseNodeData< Real, Degree > solution;

  {
  DenseNodeData< Real, Degree > constraints;
  InterpolationInfo* iInfo = NULL;
  int solveDepth = MaxSolveDepthVal;

  tree.resetNodeIndices();

  // Get the kernel density estimator [If discarding, compute anew. Otherwise, compute once.]
  {
  //profiler.start();
  density = tree.template setDensityEstimator< WEIGHT_DEGREE >(*samples, KernelDepthVal, SamplesPerNodeVal);
  //profiler.dumpOutput2(comments, "#   Got kernel density:");
  }

  // Transform the Hermite samples into a vector field [If discarding, compute anew. Otherwise, compute once.]
  {
  //profiler.start();
  normalInfo = new SparseNodeData< Point3D< Real >, NORMAL_DEGREE >();
  *normalInfo = tree.template setNormalField< NORMAL_DEGREE >(*samples, *density, pointWeightSum, BOUNDARY_NEUMANN == BOUNDARY_NEUMANN);
  //profiler.dumpOutput2(comments, "#     Got normal field:");
  }

  if (!DensityFlag) delete density, density = NULL;

  // Trim the tree and prepare for multigrid
  {
  //profiler.start();
  std::vector< int > indexMap;

  constexpr int MAX_DEGREE = NORMAL_DEGREE > Degree ? NORMAL_DEGREE : Degree;
  tree.template inalizeForBroodedMultigrid< MAX_DEGREE, Degree, BOUNDARY_NEUMANN >(FullDepthVal, typename Octree< Real >::template HasNormalDataFunctor< NORMAL_DEGREE >(*normalInfo), &indexMap);

  if (normalInfo) normalInfo->remapIndices(indexMap);
  if (density) density->remapIndices(indexMap);
  //profiler.dumpOutput2(comments, "#       Finalized tree:");
  }

  // Add the FEM constraints
  {
  //profiler.start();
  constraints = tree.template initDenseNodeData< Degree >();
  tree.template addFEMConstraints< Degree, BOUNDARY_NEUMANN, NORMAL_DEGREE, BOUNDARY_NEUMANN >(FEMVFConstraintFunctor< NORMAL_DEGREE, BOUNDARY_NEUMANN, Degree, BOUNDARY_NEUMANN >(1., 0.), *normalInfo, constraints, solveDepth);
  //profiler.dumpOutput2(comments, "#  Set FEM constraints:");
  }

  // Free up the normal info [If we don't need it for subseequent iterations.]
  delete normalInfo, normalInfo = NULL;

  // Add the interpolation constraints
  if (PointWeightVal > 0)
    {
    //profiler.start();
    iInfo = new InterpolationInfo(tree, *samples, targetValue, AdaptiveExponentVal, (Real)PointWeightVal * pointWeightSum, (Real)0);
    tree.template addInterpolationConstraints< Degree, BOUNDARY_NEUMANN >(*iInfo, constraints, solveDepth);
    //profiler.dumpOutput2(comments, "#Set point constraints:");
    }

  DumpOutput("Leaf Nodes / Active Nodes / Ghost Nodes: %d / %d / %d\n", (int)tree.leaves(), (int)tree.nodes(), (int)tree.ghostNodes());
  DumpOutput("Memory Usage: %.3f MB\n", float(MemoryInfo::Usage()) / (1 << 20));

  // Solve the linear system
  {
  //profiler.start();
  typename Octree< Real >::SolverInfo solverInfo;
  solverInfo.cgDepth = CGDepthVal, solverInfo.iters = ItersVal, solverInfo.cgAccuracy = CSSolverAccuracyVal, solverInfo.verbose = VerboseFlag, solverInfo.showResidual = ShowResidualFlag, solverInfo.lowResIterMultiplier = std::max< double >(1., LowResIterMultiplierVal);
  solution = tree.template solveSystem< Degree, BOUNDARY_NEUMANN >(FEMSystemFunctor< Degree, BOUNDARY_NEUMANN >(0, 1., 0), iInfo, constraints, solveDepth, solverInfo);
  //profiler.dumpOutput2(comments, "# Linear system solved:");
  if (iInfo) delete iInfo, iInfo = NULL;
  }
  }

  typedef PlyColorVertex< float > Vertex;

  CoredFileMeshData< Vertex > mesh;

  {
  //profiler.start();
  double valueSum = 0, weightSum = 0;
  typename Octree< Real >::template MultiThreadedEvaluator< Degree, BOUNDARY_NEUMANN > evaluator(&tree, solution, ThreadsVal);
#pragma omp parallel for num_threads( ThreadsVal ) reduction( + : valueSum , weightSum )
  for (int j = 0; j < samples->size(); j++)
    {
    ProjectiveData< OrientedPoint3D< Real >, Real >& sample = (*samples)[j].sample;
    Real w = sample.weight;
    if (w > 0) weightSum += w, valueSum += evaluator.value(sample.data.p / sample.weight, omp_get_thread_num(), (*samples)[j].node) * w;
    }
  Real isoValue = (Real)(valueSum / weightSum);
  //		if( samples ) delete samples , samples = NULL;
  //profiler.dumpOutput("Got average:");
  DumpOutput("Iso-Value: %e\n", isoValue);

  //profiler.start();
  SparseNodeData< ProjectiveData< Point3D< Real >, Real >, DATA_DEGREE >* colorData = NULL;
  if (sampleData)
    {
    colorData = new SparseNodeData< ProjectiveData< Point3D< Real >, Real >, DATA_DEGREE >();
    *colorData = tree.template setDataField< DATA_DEGREE, false >(*samples, *sampleData, (DensityEstimator*)NULL);
    delete sampleData, sampleData = NULL;
    for (const OctNode< TreeNodeData >* n = tree.tree().nextNode(); n; n = tree.tree().nextNode(n))
      {
      ProjectiveData< Point3D< Real >, Real >* clr = (*colorData)(n);
      if (clr)
        (*clr) *= (Real)pow(ColorVal, tree.depth(n));
      }
    }
  tree.template getMCIsoSurface< Degree, BOUNDARY_NEUMANN, WEIGHT_DEGREE, DATA_DEGREE >(density, colorData, solution, isoValue, mesh, !LinearFitFlag, !NonManifoldFlag, false /*PolygonMesh.set*/);
  DumpOutput("Vertices / Polygons: %d / %d\n", mesh.outOfCorePointCount() + mesh.inCorePoints.size(), mesh.polygonCount());
  //profiler.dumpOutput2(comments, "#        Got triangles:");

  if (colorData) delete colorData, colorData = NULL;
  }

  //        FreePointer( solution );


  mesh.resetIterator();
  int nr_vertices = int(mesh.outOfCorePointCount() + mesh.inCorePoints.size());
  int nr_faces = mesh.polygonCount();

  vertices.clear();
  triangles.clear();
  vertex_colors.clear();

  for (int i = 0; i<int(mesh.inCorePoints.size()); ++i)
    {
    Vertex vertex = mesh.inCorePoints[i];

    Point3D<Real> pp = iXForm * vertex.point;

    vertices.push_back(jtk::vec3<float>(pp[0], pp[1], pp[2]));

    vertex_colors.push_back(jtk::vec3<float>(vertex.color[0] / 255.f, vertex.color[1] / 255.f, vertex.color[2] / 255.f));
    }
  for (int i = 0; i < mesh.outOfCorePointCount(); ++i)
    {
    Vertex vertex;
    mesh.nextOutOfCorePoint(vertex);
    Point3D<Real> pp = iXForm * vertex.point;

    vertices.push_back(jtk::vec3<float>(pp[0], pp[1], pp[2]));
    vertex_colors.push_back(jtk::vec3<float>(vertex.color[0] / 255.f, vertex.color[1] / 255.f, vertex.color[2] / 255.f));
    }
  std::vector< CoredVertexIndex > polygon;
  while (mesh.nextPolygon(polygon))
    {
    assert(polygon.size() == 3);
    int indV[3];
    for (int i = 0; i<int(polygon.size()); i++)
      {
      if (polygon[i].inCore) indV[i] = polygon[i].idx;
      else                    indV[i] = polygon[i].idx + int(mesh.inCorePoints.size());
      }
    jtk::vec3<uint32_t> tria(indV[0], indV[1], indV[2]);
    triangles.push_back(tria);
    }

  if (density) delete density, density = NULL;
  DumpOutput2(comments, "#          Total Solve: %9.1f (s), %9.1f (MB)\n", Time() - startTime, tree.maxMemoryUsage());




  }
#pragma warning(pop)