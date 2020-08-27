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
#include <Windows.h>
#include <Psapi.h>
#include "poisson_screened/mytime.h"
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
  int MaxSolveDepthVal=-1;
  int KernelDepthVal=-1;
  int MinDepthVal=0;
  int FullDepthVal=5;
  float SamplesPerNodeVal=1.5f;
  float ScaleVal=1.1f;
  bool ConfidenceFlag=false;
  bool CleanFlag=false;
  bool DensityFlag=false;
  float PointWeightVal=4.f;
  int AdaptiveExponentVal=1;
  int BoundaryTypeVal=1;
  bool CompleteFlag=false;
  bool NonManifoldFlag=false;
  bool ShowResidualFlag=false;
  int CGDepthVal=0;
  int ItersVal=8;
  float CSSolverAccuracyVal=1e-3f;

  bool VerboseFlag=true;
  int ThreadsVal=omp_get_num_procs();
  bool LinearFitFlag=false;
  float LowResIterMultiplierVal=1.f;
  float ColorVal=16.f;

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
  //profiler.start();
  //		PointStream* pointStream;

  //		char* ext = GetFileExtension( In.value );
  //		if( Color.set && Color.value>0 )
  //		{
  //			sampleData = new std::vector< ProjectiveData< Point3D< Real > , Real > >();
  //			if     ( !strcasecmp( ext , "bnpts" ) ) pointStream = new BinaryOrientedPointStreamWithData< Real , Point3D< Real > , float , Point3D< unsigned char > >( In.value );
  //			else if( !strcasecmp( ext , "ply"   ) ) pointStream = new    PLYOrientedPointStreamWithData< Real , Point3D< Real > >( In.value , ColorInfo< Real >::PlyProperties , 6 , ColorInfo< Real >::ValidPlyProperties );
  //			else                                    pointStream = new  ASCIIOrientedPointStreamWithData< Real , Point3D< Real > >( In.value , ColorInfo< Real >::ReadASCII );
  //		}
  //		else
  //		{
  //			if     ( !strcasecmp( ext , "bnpts" ) ) pointStream = new BinaryOrientedPointStream< Real , float >( In.value );
  //			else if( !strcasecmp( ext , "ply"   ) ) pointStream = new    PLYOrientedPointStream< Real >( In.value );
  //			else                                    pointStream = new  ASCIIOrientedPointStream< Real >( In.value );
  //		}
  //		delete[] ext;
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
  /*
  for (int i = 0; i < nr_faces; ++i)
    {
    PlyFace ply_face;
    mesh.nextPolygon(polygon);
    ply_face.nr_vertices = int(polygon.size());
    ply_face.vertices = new int[polygon.size()];
    for (int i = 0; i<int(polygon.size()); ++i)
      {
      if (polygon[i].inCore)
        ply_face.vertices[i] = polygon[i].idx;
      else
        ply_face.vertices[i] = polygon[i].idx + int(mesh.inCorePoints.size());
      }
    for (int j = 2; j < ply_face.nr_vertices; ++j)
      {
      jtk::vec3<uint32_t> tria(ply_face.vertices[0], ply_face.vertices[1], ply_face.vertices[j]);
      triangles.push_back(tria);
      }
    delete[] ply_face.vertices;
    }
    */
  /*
  int vm = mesh.outOfCorePointCount() + mesh.inCorePoints.size();
  for (auto pt = mesh.inCorePoints.begin(); pt != mesh.inCorePoints.end(); ++pt)
    {
    Point3D<Real> pp = iXForm * pt->point;
    vcg::tri::Allocator<CMeshO>::AddVertex(pm, Point3m(pp[0], pp[1], pp[2]));
    pm.vert.back().Q() = pt->value;
    pm.vert.back().C()[0] = pt->color[0];
    pm.vert.back().C()[1] = pt->color[1];
    pm.vert.back().C()[2] = pt->color[2];
    }
  for (int ii = 0; ii < mesh.outOfCorePointCount(); ii++) {
    Vertex pt;
    mesh.nextOutOfCorePoint(pt);
    Point3D<Real> pp = iXForm * pt.point;
    vcg::tri::Allocator<CMeshO>::AddVertex(pm, Point3m(pp[0], pp[1], pp[2]));
    pm.vert.back().Q() = pt.value;
    pm.vert.back().C()[0] = pt.color[0];
    pm.vert.back().C()[1] = pt.color[1];
    pm.vert.back().C()[2] = pt.color[2];
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
    vcg::tri::Allocator<CMeshO>::AddFace(pm, &pm.vert[indV[0]], &pm.vert[indV[1]], &pm.vert[indV[2]]);
    }
 
 */
  //		if( colorData ) delete colorData , colorData = NULL;


  if (density) delete density, density = NULL;
  DumpOutput2(comments, "#          Total Solve: %9.1f (s), %9.1f (MB)\n", Time() - startTime, tree.maxMemoryUsage());



  /*
  int MinDepth = 0;
  int Depth = par.depth;
  int FullDepth = par.octree_depth;
  int kernelDepth = par.depth - 2;
  double SamplesPerNode = double(par.samples_per_node);
  double Scale = 1.1;
  bool Confidence = false;
  bool NormalWeights = false;
  double PointWeight = 4.0;
  int AdaptiveExponent = 1;
  int BoundaryType = 1;
  bool Complete = false;
  auto xForm = XForm4x4< double >::Identity();
  bool Density = false;
  int MaxSolveDepth = Depth;
  bool ShowResidual = false;
  int Iters = 8;
  int CGDepth = par.solver_divide;
  double CSSolverAccuracy = 1e-3;

  Reset<double>();
  Octree< double > tree;
  tree.threads = omp_get_num_procs();
  OctNode< TreeNodeData >::SetAllocator(MEMORY_ALLOCATOR_BLOCK_SIZE);

  double maxMemoryUsage;
  double t = Time();
  double tt = Time();
  tree.maxMemoryUsage = 0;

  Octree< double >::PointInfo* pointInfo = new Octree< double >::PointInfo();
  Octree< double >::NormalInfo* normalInfo = new Octree< double >::NormalInfo();
  std::vector< double >* kernelDensityWeights = new std::vector< double >();
  std::vector< double >* centerWeights = new std::vector< double >();
  PointStream< double >* pointStream = new MyPointStream(pts, normals);
  pointStream->reset();
  int pointCount = tree.template SetTree< double >(pointStream, MinDepth, Depth, FullDepth, kernelDepth, SamplesPerNode, Scale, Confidence, NormalWeights, PointWeight, AdaptiveExponent, *pointInfo, *normalInfo, *kernelDensityWeights, *centerWeights, BoundaryType, xForm, Complete);

  if (!Density)
    {
    delete kernelDensityWeights;
    kernelDensityWeights = nullptr;
    }

  maxMemoryUsage = tree.maxMemoryUsage;
  t = Time();
  tree.maxMemoryUsage = 0;
  Pointer(double) constraints = tree.SetLaplacianConstraints(*normalInfo);
  delete normalInfo;
  maxMemoryUsage = std::max< double >(maxMemoryUsage, tree.maxMemoryUsage);

  t = Time();
  tree.maxMemoryUsage = 0;
  Pointer(double) solution = tree.SolveSystem(*pointInfo, constraints, ShowResidual, Iters, MaxSolveDepth, CGDepth, CSSolverAccuracy);
  delete pointInfo;
  FreePointer(constraints);
  maxMemoryUsage = std::max< double >(maxMemoryUsage, tree.maxMemoryUsage);

  CoredFileMeshData< PlyVertex< float > > mesh;
  tree.maxMemoryUsage = 0;
  //float isoValue = tree.GetIsoValue(solution, *centerWeights);
  double isoValue = 0.0;
  delete centerWeights;
  tree.GetMCIsoSurface(kernelDensityWeights ? GetPointer(*kernelDensityWeights) : NullPointer< double >(), solution, isoValue, mesh, true, true, false);

  auto iXForm = xForm.inverse();

  //PlyWritePolygons("PoissonMesh.ply", &mesh, PLY_ASCII, NULL, 0, iXForm);

  mesh.resetIterator();

  int nr_vertices = int(mesh.outOfCorePointCount() + mesh.inCorePoints.size());
  int nr_faces = mesh.polygonCount();

  vertices.clear();
  triangles.clear();

  for (int i = 0; i<int(mesh.inCorePoints.size()); ++i)
    {
    PlyVertex< float > vertex = mesh.inCorePoints[i];
    vertices.push_back(jtk::vec3<float>(vertex.point[0], vertex.point[1], vertex.point[2]));
    }
  for (int i = 0; i < mesh.outOfCorePointCount(); ++i)
    {
    PlyVertex< float > vertex;
    mesh.nextOutOfCorePoint(vertex);
    vertices.push_back(jtk::vec3<float>(vertex.point[0], vertex.point[1], vertex.point[2]));
    }
  std::vector< CoredVertexIndex > polygon;
  for (int i = 0; i < nr_faces; ++i)
    {
    PlyFace ply_face;
    mesh.nextPolygon(polygon);
    ply_face.nr_vertices = int(polygon.size());
    ply_face.vertices = new int[polygon.size()];
    for (int i = 0; i<int(polygon.size()); ++i)
      {
      if (polygon[i].inCore)
        ply_face.vertices[i] = polygon[i].idx;
      else
        ply_face.vertices[i] = polygon[i].idx + int(mesh.inCorePoints.size());
      }
    for (int j = 2; j < ply_face.nr_vertices; ++j)
      {
      jtk::vec3<uint32_t> tria(ply_face.vertices[0], ply_face.vertices[1], ply_face.vertices[j]);
      triangles.push_back(tria);
      }
    delete[] ply_face.vertices;
    }
    */
  }

#pragma warning(pop)