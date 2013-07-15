#ifndef CUBICBEZIERSURFACE_HH
#define CUBICBEZIERSURFACE_HH

#include <math/Matrix4x4.hh>
#include <math/Vector3D.hh>
#include <geometry/CubicBezierCurve3D.hh>

namespace geometry{

class CubicBezierSurface
{
public:
	/*! The control points should be arranged as: 4 points for 1st row, 4 points for 2nd row...
	*/
	CubicBezierSurface(Vector3D controlPoints[16]);

	CubicBezierSurface(CubicBezierSurface & bezierCurve);

	void Copy(CubicBezierSurface & bezierCurve);

	CubicBezierSurface operator =(CubicBezierSurface & bezierCurve){ Copy(bezierCurve); return *this; };

	Vector3D GetPoint(double u, double v);

	/*! Reverse the order of control points
	\param direction u(=0) or v(=1)
	*/
	void Reverse(int direction);

	/*! Get the 16 control points in array form.
	*/
	Vector3D * GetControlPoints(){ return _controlPoints; };

    void SetPSpace( const PSpace& nPSpace ) { _pSpace = nPSpace; };
	PSpace & GetPSpace(){return _pSpace;};

	/*! Get one of the 4 boundary curves
	\param boundary index (0 to 3) 0:(u=0.0), 1:(u=1.0), 2:(v=0.0), 3(v=1.0)
	*/
	CubicBezierCurve3D GetBoundaryCurve(int index);

private:
	Vector3D _controlPoints[16];
	PSpace _pSpace;
};

class QuadTree
{
public:
	QuadTree(QuadTree * parent);

	~QuadTree(){
		Clean();
	}

	void Clean();

	QuadTree * GetParent(){ return _parent; };

	QuadTree ** GetChildren(){ return _children; };

	QuadTree ** GetNeighbors(){ return _neighbors; };

	int * GetRowRange(){ return _rowRange; };

	int * GetColRange(){ return _colRange; };

	virtual void MakeChildren();

	bool HaveChildren(){ return _children[0] == NULL ? false : true; };

	bool HaveNeighbors(){ return _neighbors[0] == NULL ? false : true; };

	bool IsRoot(){return _parent == NULL ? true : false; };

	int GetPosition();

	QuadTree * GetRoot(){ return _root; };

	class QuadTreeIterator
	{		
	public:
		QuadTreeIterator(QuadTree * quadTree, int depth = 0);

		QuadTree * operator *() const;

		bool operator == (const QuadTreeIterator & iter) const{
			return _current == (*iter) ? true : false;
		}

		bool operator != (const QuadTreeIterator & iter) const{
			return !(*this == iter);
		}

		bool operator < (const QuadTreeIterator & iter) const{
			return _current < (*iter);
		}

		QuadTreeIterator operator = (QuadTree * quadTree){
			_current = quadTree;
			_position = quadTree->GetPosition();
			_depth = 0;
			return *this;
		};

		QuadTreeIterator operator ++();

		QuadTreeIterator Next();

		QuadTreeIterator Up();

		QuadTreeIterator GoTo(int row, int col, int maxDepth = -1);

		int GetDepth(){ return _depth; };

	private:
		QuadTree * _current;

		int _position;
		int _depth;
	};

private:
	QuadTree * _parent;
	QuadTree * _root;

	//! Child branches
	// 2 3
	// 0 1
	QuadTree * _children[4];

	//! Neighbor branches
	// * 3 *
	// 0 C 1
	// * 2 *
	QuadTree * _neighbors[4];

	//! Note: First and last row and column ranges overlap with the neighbors!
	int _rowRange[2];

	//! Note: First and last row and column ranges overlap with the neighbors!
	int _colRange[2];
};

class BezierQuadTree : public QuadTree
{
public:
	BezierQuadTree();

	BezierQuadTree(BezierQuadTree * parent);

	~BezierQuadTree(){
		Clean();
	}

	void Clean();

	void MakeChildren();

	CubicBezierSurface * GetSurface(){ return _haveSurface ? _surface : NULL; };

	void SetSurface(CubicBezierSurface * surface);

	double GetVariance(){ return _variance; };

	void SetVariance(double variance){ _variance = variance; };

	bool HaveSurface(){ return _haveSurface; };

	double GridToParameters(const double x, const double y, int index[2], double param[2], int maxDepth = -1);

private:
	CubicBezierSurface * _surface;
	double _variance;
	bool _haveSurface;

	double Brent(CubicBezierSurface * surface, int i, double x, double y, double x1, double x2, int maxIteration, double tolerance);
};

class CubicBezierSurfaceUtils
{
public:
	static CubicBezierSurface Divide(CubicBezierSurface & bezierCurve, double uRange[2], double vRange[2]);

	static CubicBezierSurface FitSurfaceToGrid(Grid2D & grid, int rowRange[2], int colRange[2], int sampleInterval[2]);

	static void FitBoundaryToGrid(CubicBezierSurface & bezier, Grid2D & grid, int dir, int n, int rowRange[2], int colRange[2], int sampleInterval);

	static void AdjustInternalControlPointsToFitGrid(CubicBezierSurface & bezier, Grid2D & grid, int rowRange[2], int colRange[2], int sampleInterval[2]);

	/*!
	\param dir direction (0 = row, 1 = col)
	\param n = n-th row or column
	*/
	static void ExtractGridPolyline(Grid2D & grid, int dir, int n, int range[2], int interval, std::vector<Vector3D> & polyline);

	static double GetHeightVariance(Grid2D & grid, int rowRange[2], int colRange[2], int sampleInterval[2]);

	static void DumpInventorFile(std::ostream& ovstr, CubicBezierSurface & bezier, float color[3]);

	static void DumpPoints(std::ostream& ovstr, CubicBezierSurface & bezierCurve, double uRange[2], double vRange[2], double uInterval, double vInterval);

	static void DumpInventorPoints(std::ostream& ovstr, CubicBezierSurface & bezierCurve, double uRange[2], double vRange[2], double uInterval, double vInterval, float color[3], double pointSize = 1.0);

	static void DumpInventorControlPoints(std::ostream& ovstr, CubicBezierSurface & bezierSurface, float color[3], double pointSize);

	static void DumpInventorFile(std::ostream& ovstr, BezierQuadTree * root, float color[3], int depthLimit = -1);

	static void DumpControlPointsInventorFile(std::ostream& ovstr, BezierQuadTree * root, float color[3], float pointSize = 1.0, int depthLimit = -1);

	static void DumpPatchBoundariesInventorFile(std::ostream& ovstr, BezierQuadTree * root, float color[3], float lineThickness = 1.0, int depthLimit = -1);

	static BezierQuadTree * BuildBezierPatchesFromGrid(Grid2D & grid, double resolution, bool doSmoothing = true, bool doLocalFitting = true);
private:
	static BezierQuadTree * makeQuadTree(Grid2D & grid, double resolution, int depthRange[2]);
	static void fitPatchToGrid(BezierQuadTree::QuadTreeIterator & iter, Grid2D & grid, std::set< int > & gridIndexToSmooth, int shallowestDepth);
	static void smoothChildPatches(BezierQuadTree::QuadTreeIterator & iter, Grid2D & grid);
	static void doGlobalSmoothing(BezierQuadTree * root, Grid2D &grid, std::set< int > & gridIndexToSmooth, int depthLimit = -1);
	static void smoothPoints(std::map<int, std::vector< std::vector < Vector3D * > > > & nodeToControlPoints);
	static void alignAllPatchesToGrid(BezierQuadTree * root, int depth = -1);
	static void alignPatchToGrid(BezierQuadTree * patch);
	static void makeChildSubdivisionSurface(BezierQuadTree * bTree);
	static std::set< int >getDepthsAroundGrid(Grid2D &grid, BezierQuadTree * root, int row, int col, int upperLimit = -1, int lowerLimit = -1);
	static std::set< CubicBezierSurface * >getBezierSurfacesAroundGrid(Grid2D &grid, BezierQuadTree * root, int row, int col);
	static Matrix4X4 invertMatrix4X4(Matrix4X4 & m);
};

} // namespace geometry

#endif
