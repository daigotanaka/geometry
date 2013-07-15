#include "CubicBezierSurface.hh"

#include <iostream>
#include <fstream>

using namespace std;

namespace geometry{

CubicBezierSurface::CubicBezierSurface(Vector3D controlPoints[16])
{
	for (int i = 0; i < 16; i++){
		_controlPoints[i] = controlPoints[i];
	}
}

CubicBezierSurface::CubicBezierSurface(CubicBezierSurface & bezier)
{
	Copy(bezier);
}

void
CubicBezierSurface::Copy(CubicBezierSurface & bezier)
{
	for (int i = 0; i < 16; i++){
		_controlPoints[i] = bezier.GetControlPoints()[i];
	}
	SetPSpace( bezier.GetPSpace() );
}

Vector3D
CubicBezierSurface::GetPoint(double u, double v)
{
	Vector3D point(0, 0, 0);
	double bi[4] = {
		CubicBezierCurve3D::Bernstein(u, 0),
		CubicBezierCurve3D::Bernstein(u, 1),
		CubicBezierCurve3D::Bernstein(u, 2),
		CubicBezierCurve3D::Bernstein(u, 3),
	};
	double bj[4] = {
		CubicBezierCurve3D::Bernstein(v, 0),
		CubicBezierCurve3D::Bernstein(v, 1),
		CubicBezierCurve3D::Bernstein(v, 2),
		CubicBezierCurve3D::Bernstein(v, 3),
	};
	for (int j = 0; j < 4; j++){
		Vector3D b0 = bj[j] * bi[0] * _controlPoints[j * 4 + 0];
		Vector3D b1 = bj[j] * bi[1] * _controlPoints[j * 4 + 1];
		Vector3D b2 = bj[j] * bi[2] *_controlPoints[j * 4 + 2];
		Vector3D b3 = bj[j] * bi[3] * _controlPoints[j * 4 + 3];
		point = point + b0 + b1 + b2 + b3;
	}

	return point;
}

void
CubicBezierSurface::Reverse(int direction)
{
	swap(_controlPoints[0], _controlPoints[3]);
	swap(_controlPoints[1], _controlPoints[2]);
}

CubicBezierCurve3D
CubicBezierSurface::GetBoundaryCurve(int index)
{
	Vector3D controlPoints[4];
	if (index < 2){
		for (int i = 0; i < 4; i++){
			controlPoints[i] = this->GetControlPoints()[index * 12 + i];
		}
	}
	else{
		index -= 2;
		for (int i = 0; i < 4; i++){
			controlPoints[i] = this->GetControlPoints()[i * 4 + index * 3];
		}
	}
	return CubicBezierCurve3D(controlPoints);
}

QuadTree::QuadTree(QuadTree * parent):
_parent(parent)
{
	_children[0] = NULL;
	_neighbors[0] = NULL;
	if (parent != NULL){
		_root = _parent->GetRoot();
	}
	else{
		_root = this;
	}
}

void
QuadTree::Clean()
{
	if (_children[0] == NULL){
		return;
	}

	for (int i = 0; i < 4; i++){
		delete _children[i];
	}
	_children[0] = NULL;
	_neighbors[0] = NULL;
}

void
QuadTree::MakeChildren()
{
	if (HaveChildren()){
		return;
	}

	if (_rowRange[1] - _rowRange[0] < 2 || _colRange[1] - _colRange[0] < 2){
		return;
	}

	int * rowRange;
	int * colRange;

	_children[0] = new QuadTree(this);
	rowRange = _children[0]->GetRowRange();
	colRange = _children[0]->GetColRange();
	rowRange[0] = _rowRange[0];
	rowRange[1] = (int)(_rowRange[0] + 0.5 * (_rowRange[1] - _rowRange[0]));
	colRange[0] = _colRange[0];
	colRange[1] = (int)(_colRange[0] + 0.5 * (_colRange[1] - _colRange[0]));

	_children[1] = new QuadTree(this);
	rowRange = _children[1]->GetRowRange();
	colRange = _children[1]->GetColRange();
	rowRange[0] = _rowRange[0];
	rowRange[1] = (int)(_rowRange[0] + 0.5 * (_rowRange[1] - _rowRange[0]));
	colRange[0] = (int)(_colRange[0] + 0.5 * (_colRange[1] - _colRange[0]));
	colRange[1] = _colRange[1];

	_children[2] = new QuadTree(this);
	rowRange = _children[2]->GetRowRange();
	colRange = _children[2]->GetColRange();
	rowRange[0] = (int)(_rowRange[0] + 0.5 * (_rowRange[1] - _rowRange[0]));
	rowRange[1] = _rowRange[1];
	colRange[0] = _colRange[0];
	colRange[1] = (int)(_colRange[0] + 0.5 * (_colRange[1] - _colRange[0]));

	_children[3] = new QuadTree(this);
	rowRange = _children[3]->GetRowRange();
	colRange = _children[3]->GetColRange();
	rowRange[0] = (int)(_rowRange[0] + 0.5 * (_rowRange[1] - _rowRange[0]));
	rowRange[1] = _rowRange[1];
	colRange[0] = (int)(_colRange[0] + 0.5 * (_colRange[1] - _colRange[0]));
	colRange[1] = _colRange[1];
}

int
QuadTree::GetPosition()
{
	if (IsRoot()){
		return -1;
	}

	QuadTree ** siblings = GetParent()->GetChildren();

	for (int i = 0; i < 4; i++){
		if (siblings[i] == this){
			return i;
		}
	}
	return -1;
}

QuadTree::QuadTreeIterator::QuadTreeIterator(QuadTree * quadTree, int depth):
_current(quadTree), _depth(depth)
{
	if (quadTree != NULL){
		_position = _current->GetPosition();
	}
}

QuadTree * 
QuadTree::QuadTreeIterator::operator *() const
{
	return _current;
}

QuadTree::QuadTreeIterator
QuadTree::QuadTreeIterator::operator ++()
{
	// go deeper if we can
	if (_current->HaveChildren()){
		_current = _current->GetChildren()[0];
		_position = 0;
		_depth++;

		return *this;
	}

	return Next();
}

QuadTree::QuadTreeIterator
QuadTree::QuadTreeIterator::Next()
{
	if (_current == NULL){
		return NULL;
	}

	QuadTree * parent = _current->GetParent();

	// try to find the next siblings
	_position++;
	if (0 < _position && _position < 4){
		_current = parent->GetChildren()[_position];
		return *this;
	}

	// no more slblings. Move to parents' sibling
	if (_position == 4){
		Up();
		return Next();
	}

	_current = NULL;
	_position = -1;
	return _current;
}

QuadTree::QuadTreeIterator
QuadTree::QuadTreeIterator::Up()
{
	_current = _current->GetParent();
	_position = _current->GetPosition();
	_depth--;
	return _current;
}

QuadTree::QuadTreeIterator
QuadTree::QuadTreeIterator::GoTo(int row, int col, int maxDepth)
{
	int * rowRange = _current->GetRowRange();
	int * colRange = _current->GetColRange();
	if (row < rowRange[0]  || rowRange[1] < row
		|| col < colRange[0] || colRange[1] < col){
			_current = NULL;
			_position = -1;
			_depth = -1;
			return NULL;
	}

	while(_current->HaveChildren() && ( maxDepth == -1 || _depth < maxDepth) ){
		rowRange = _current->GetRowRange();
		colRange = _current->GetColRange();
		_depth++;
		if (row <= (int)(rowRange[0] + 0.5 * (rowRange[1] - rowRange[0])) && col <= (int)(colRange[0] + 0.5 * (colRange[1] - colRange[0]))){
			_current = _current->GetChildren()[0];
		}
		else if (row <= (int)(rowRange[0] + 0.5 * (rowRange[1] - rowRange[0]))){
			_current = _current->GetChildren()[1];
		}
		else if (col <= (int)(colRange[0] + 0.5 * (colRange[1] - colRange[0]))){
			_current = _current->GetChildren()[2];
		}
		else{
			_current = _current->GetChildren()[3];
		}
	}
	_position = _current->GetPosition();

	return *this;
}

BezierQuadTree::BezierQuadTree():
QuadTree(NULL), _surface(NULL), _haveSurface(false)
{
};

BezierQuadTree::BezierQuadTree(BezierQuadTree * parent):
QuadTree(parent), _surface(NULL), _haveSurface(false)
{
};

void
BezierQuadTree::MakeChildren()
{
	if (HaveChildren()){
		return;
	}

	int * parentRowRange = GetRowRange();
	int * parentColRange = GetColRange();

	if (parentRowRange[1] - parentRowRange[0] < 2 || parentColRange[1] - parentColRange[0] < 2){
		return;
	}

	int * rowRange;
	int * colRange;

	GetChildren()[0] = new BezierQuadTree(this);
	rowRange = GetChildren()[0]->GetRowRange();
	colRange = GetChildren()[0]->GetColRange();
	rowRange[0] = parentRowRange[0];
	rowRange[1] = (int)(parentRowRange[0] + 0.5 * (parentRowRange[1] - parentRowRange[0]));
	colRange[0] = parentColRange[0];
	colRange[1] = (int)(parentColRange[0] + 0.5 * (parentColRange[1] - parentColRange[0]));

	GetChildren()[1] = new BezierQuadTree(this);
	rowRange = GetChildren()[1]->GetRowRange();
	colRange = GetChildren()[1]->GetColRange();
	rowRange[0] = parentRowRange[0];
	rowRange[1] = (int)(parentRowRange[0] + 0.5 * (parentRowRange[1] - parentRowRange[0]));
	colRange[0] = (int)(parentColRange[0] + 0.5 * (parentColRange[1] - parentColRange[0]));
	colRange[1] = parentColRange[1];

	GetChildren()[2] = new BezierQuadTree(this);
	rowRange = GetChildren()[2]->GetRowRange();
	colRange = GetChildren()[2]->GetColRange();
	rowRange[0] = (int)(parentRowRange[0] + 0.5 * (parentRowRange[1] - parentRowRange[0]));
	rowRange[1] = parentRowRange[1];
	colRange[0] = parentColRange[0];
	colRange[1] = (int)(parentColRange[0] + 0.5 * (parentColRange[1] - parentColRange[0]));

	GetChildren()[3] = new BezierQuadTree(this);
	rowRange = GetChildren()[3]->GetRowRange();
	colRange = GetChildren()[3]->GetColRange();
	rowRange[0] = (int)(parentRowRange[0] + 0.5 * (parentRowRange[1] - parentRowRange[0]));
	rowRange[1] = parentRowRange[1];
	colRange[0] = (int)(parentColRange[0] + 0.5 * (parentColRange[1] - parentColRange[0]));
	colRange[1] = parentColRange[1];
}

void
BezierQuadTree::Clean()
{
	delete _surface;
	_haveSurface = false;
}

void
BezierQuadTree::SetSurface(CubicBezierSurface * surface)
{
	if  (surface != NULL) {
		_surface = surface;
		_haveSurface = true;
	}

	Vector3D * controlPoints = _surface->GetControlPoints();
	Vector3D(controlPoints[0][0], controlPoints[0][1], controlPoints[0][2]);
}

double
BezierQuadTree::GridToParameters(const double row, const double col, int index[2], double param[2], int maxDepth)
{
	double tolerance = TOLERANCE;

	int * rowRange = GetRowRange();
	int * colRange = GetColRange();

	if (col < (double)colRange[0] || (double) colRange[1] < col
		|| row < (double)rowRange[0] || (double) rowRange[1] < row
		){
			return (double)std::numeric_limits<double>::max();
	}

	index[0] = (int)floor(row); // row
	index[1] = (int)floor(col); // col
	double x = col * this->GetSurface()->GetPSpace().GetAxisLength( 0 );
	double y = row * this->GetSurface()->GetPSpace().GetAxisLength( 1 );
	Vector3D gridPoint(x, y, 0.0);

	BezierQuadTree::QuadTreeIterator iter((QuadTree *)this);
	iter.GoTo(index[0], index[1], maxDepth);
	BezierQuadTree * bTree = (BezierQuadTree *)(*iter);
	rowRange = bTree->GetRowRange();
	colRange = bTree->GetColRange();
	CubicBezierSurface * surface = bTree->GetSurface();

	param[0] = (col - colRange[0]) / (colRange[1] - colRange[0]); // u
	param[1] = (row - rowRange[0]) / (rowRange[1] - rowRange[0]); // v
	for (int i = 0; i < 2; i++){
		if (param[i] < 0.0){
			param[i] = 0.0;
		}
		if (param[i] > 1.0){
			param[i] = 1.0;
		}
	}

	double error[2];
	Vector3D point;

	for (int iteration = 0; iteration < 2; iteration++){
		point = surface->GetPoint(param[0], param[1]);
		for (int i = 0; i < 2; i++){
			if (param[i] == 0.0 || param[i] == 1.0){
				continue;
			}
			error[i] = gridPoint[i] - point[i];

			if (abs(error[i]) > tolerance && 0.0 < param[i] && param[i] < 1.0){
				if (error[i] < 0){
					param[i] = Brent(surface, i, gridPoint[i], param[1 - i], 0.0, param[i], 20, tolerance);
				}
				else{
					param[i] = Brent(surface, i, gridPoint[i], param[1 - i], param[i], 1.0, 20, tolerance);
				}

				if (param[i] < 0){
					param[i] = 0.0;
				}
				else if (param[i] > 1.0){
					param[i] = 1.0;
				}
				point = surface->GetPoint(param[0], param[1]);
				error[0] = gridPoint[0] - point[0];
				error[1] = gridPoint[1] - point[1];
			}
		}
	}

	rowRange = GetRowRange();
	colRange = GetColRange();

	return point[2];

	/////////////////////////////
	if ( (0.0 <= param[0] && param[0] <= 1.0 && 0.0 <= param[1] && param[1] <= 1.0 )
		|| index[0] == rowRange[0]
		|| index[0] == rowRange[1]
		|| index[1] == colRange[0]
		|| index[1] == colRange[1]
		){
		return point[2];
	}

	return 10;

	bool check[9];
	for (int i = 0; i < 9; i++){
		check[i] = true;
	}
	check[4] = false;
	if (param[0] > 0.0){
		for (int i = 0; i < 3; i++){
			check[i * 3] = false;
		}
	}
	if (param[0] < 1.0){
		for (int i = 0; i < 3; i++){
			check[i * 3 + 2] = false;
		}
	}
	if (param[1] > 0.0){
		for (int i = 0; i < 3; i++){
			check[i] = false;
		}
	}
	if (param[1] < 1.0){
		for (int i = 0; i < 3; i++){
			check[6 + i] = false;
		}
	}

	int bestIndex[2] = {index[0], index[1]};
	double bestParam[2] = {param[0], param[1]};

	for (int i = -1; i <= 1; i++){ // row
		for (int j = -1; j <= 1; j++){ // col
			if (check[ (i + 1) * 3 + j + 1 ] == false){
				continue;
			}

			int currentIndex[2] = { index[0] + i, index[1] + j};

			BezierQuadTree::QuadTreeIterator iter((QuadTree *)this);
			iter.GoTo(currentIndex[0], currentIndex[1], maxDepth);
			if (iter == NULL){
				continue;
			}
			BezierQuadTree * bTree = (BezierQuadTree *)(*iter);
			CubicBezierSurface * surface = bTree->GetSurface();

			double currentParam[2];
			currentParam[0] = j == 0 ? param[0] : (1.0 - (double)j) * 0.5;
			currentParam[1] = i == 0 ? param[1] : (1.0 - (double)i) * 0.5;
			double currentError[2];

			Vector3D currentPoint;
			
			for (int iteration = 0; iteration < 2; iteration++){
			currentPoint = surface->GetPoint(currentParam[0], currentParam[1]);
			for (int k = 0; k < 2; k++){
				currentError[k] = gridPoint[k] - currentPoint[k];
				if (currentError[k] < 0){
					currentParam[k] = Brent(surface, k, gridPoint[k], currentParam[1 - k], 0.0, currentParam[k], 20, tolerance);
				}
				else{
					currentParam[k] = Brent(surface, k, gridPoint[k], currentParam[1 - k], currentParam[k], 1.0, 20, tolerance);
				}

				if (currentParam[k] < 0){
					currentParam[k] = 0.0;
				}
				else if (currentParam[k] > 1.0){
					currentParam[k] = 1.0;
				}
				currentPoint = surface->GetPoint(currentParam[0], currentParam[1]);
				currentError[0] = gridPoint[0] - currentPoint[0];
				currentError[1] = gridPoint[1] - currentPoint[1];
			}
			if (currentError[0] * currentError[0] + currentError[1] * currentError[1] < error[0] * error[0] + error[1] * error[1]){
				bestIndex[0] = currentIndex[0];
				bestIndex[1] = currentIndex[1];
				bestParam[0] = currentParam[0];
				bestParam[1] = currentParam[1];
				point = currentPoint;
			}
			}

		}
	}

	index[0] = bestIndex[0];
	index[1] = bestIndex[1];
	param[0] = bestParam[0];
	param[1] = bestParam[1];

	return point[2];
}

#define SIGN(a, b) ((b) >= 0.0 ? abs(a) : -abs(a))
double
BezierQuadTree::Brent(CubicBezierSurface * surface, int i, double x, double y, double x1, double x2, int maxIteration, double tol)
{
	int iter;
	double a = x1;
	double b = x2;
	double c = x2;
	double d;
	double e;
	double min1;
	double min2;

	Vector3D point;

	double param[2];
	param[1 - i] = y;

	param[i] = a;
	point = surface->GetPoint(param[0], param[1]);
	double fa = point[i] - x;

	param[i] = b;
	point = surface->GetPoint(param[0], param[1]);
	double fb = point[i] - x;

	double fc;
	double p;
	double q;
	double r;
	double s;
	double tol1;
	double xm;

	double failSafeVal = abs(fa) < abs(fb) ? x1 : x2;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)){
		thow Excetion( "BezierQuadTree::Brent Root must be bracketed in zBrent" );
	}

	fc=fb;
	for (iter=1;iter<= maxIteration;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a; // Rename a, b, c and adjust bounding interval
			fc=fa; // d.
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1 = 2.0 * std::numeric_limits<double>::epsilon() * fabs(b) + 0.5 * tol; // Convergence check.
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa; // Attempt inverse quadratic interpolation.
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q; // Check whether in bounds.
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d; // Accept interpolation.
				d=p/q;
			} else {
				d=xm; // Interpolation failed, use bisection.
				e=d;
			}
		} else { // Bounds decreasing too slowly, use bisection.
			d=xm;
			e=d;
		}
		a=b; // Move last best guess to a.
		fa=fb;
		if (fabs(d) > tol1) // Evaluate new trial root.
			b += d;
		else
			b += SIGN(tol1,xm);

		param[i] = b;
		point = surface->GetPoint(param[0], param[1]);
		fb = point[i] - x;
	}

	throw Exception( "BezierQuadTree::Brent Maximum number of iterations exceeded in zBrent" );
}

CubicBezierSurface
CubicBezierSurfaceUtils::Divide(CubicBezierSurface & bezierSurface, double uRange[2], double vRange[2])
{
	if (uRange[0] < 0.0 || uRange[1] > 1.0 || uRange[0] > uRange[1]){
		return NULL;
	}
	if (vRange[0] < 0.0 || vRange[1] > 1.0 || vRange[0] > vRange[1]){
		return NULL;
	}

	Vector3D controlPoints[16];
	for (int i = 0; i < 16; i++){
		controlPoints[i] = bezierSurface.GetControlPoints()[i];
	}

	// first divide by u direction
	if (uRange[0] != 0.0 || uRange[1] != 1.0){
		for (int i = 0; i < 4; i++){
			CubicBezierCurve3D curve(&(controlPoints[i * 4]));
			CubicBezierCurve3D newCurve = CubicBezierCurve3DUtils::Divide(curve, uRange[0], uRange[1]);
			for (int j = 0; j < 4; j++){
				controlPoints[i * 4 + j] = newCurve.GetControlPoints()[j];
			}
		}
	}

	// now divide by v direction
	if (vRange[0] != 0.0 || vRange[1] != 1.0){
		for (int i = 0; i < 4; i++){
			Vector3D transposedControlPoints[4] = {
				controlPoints[0 + i],
				controlPoints[4 + i],
				controlPoints[8 + i],
				controlPoints[12 + i]
			};
			CubicBezierCurve3D curve(transposedControlPoints);
			CubicBezierCurve3D newCurve = CubicBezierCurve3DUtils::Divide(curve, vRange[0], vRange[1]);
			for (int j = 0; j < 4; j++){
				controlPoints[j * 4 + i] = newCurve.GetControlPoints()[j];
			}
		}
	}

	CubicBezierSurface division(controlPoints);
	division.SetPSpace(bezierSurface.GetPSpace());

	return division;
}

BezierQuadTree *
CubicBezierSurfaceUtils::BuildBezierPatchesFromGrid(Grid2D & grid, double resolution, bool doSmoothing, bool doLocalFitting)
{
	int depthRange[2];
	BezierQuadTree * root = makeQuadTree(grid, resolution, depthRange);

	// Fit local patches
	// Fit the patches of the same level first and secure the continuity among siblings before going deeper.
	Vector3D cp[16];
	CubicBezierSurface * bezierSurface = new CubicBezierSurface( cp );
	bezierSurface->SetPSpace( grid.GetPSpace() );
	root->SetSurface( bezierSurface );
	BezierQuadTree::QuadTreeIterator iter = ((QuadTree *)root);

	// Fit the root surface
	set<int> gridIndexToSmooth;
	fitPatchToGrid(iter, grid, gridIndexToSmooth, depthRange[0]);

	gridIndexToSmooth.clear();

	int maxDepth = 0;
	int currentDepth = 0;
	while (currentDepth <= maxDepth){
		maxDepth = maxDepth < iter.GetDepth() ? iter.GetDepth() : maxDepth;

		if (iter.GetDepth() != currentDepth){
			if ((*iter) != NULL){
				++iter;
			}
			if ((*iter) == NULL){
				currentDepth++;
				if (doSmoothing){
					if (doLocalFitting && doSmoothing){
						doGlobalSmoothing(root, grid, gridIndexToSmooth, currentDepth);
					}
					gridIndexToSmooth.clear();
				}

				iter = root;
				++iter;
			}
			continue;
		}

		if (!(*iter)->HaveChildren()){
			++iter;
			continue;
		}

		BezierQuadTree * bTree = (BezierQuadTree *)*iter;
		makeChildSubdivisionSurface(bTree);

		// Localy fit the child patches
		if (doLocalFitting){
			BezierQuadTree::QuadTreeIterator childIter = iter;
			++childIter; // moved to the first child
			for (int childIndex = 0; childIndex < 4; childIndex++){
				fitPatchToGrid(childIter, grid, gridIndexToSmooth, depthRange[0]);
				childIter.Next();
			}
		}

		// Move on
		++iter;
	}

	return root;
}

BezierQuadTree *
CubicBezierSurfaceUtils::makeQuadTree(Grid2D & grid, double resolution, int depthRange[2])
{
	int minGridPoints = 8;
	BezierQuadTree * root = new BezierQuadTree();
	int * rowRange = root->GetRowRange();
	int * colRange = root->GetColRange();
	rowRange[0] = 0;
	rowRange[1] = grid.numRows() - 1;
	colRange[0] = 0;
	colRange[1] = grid.numColumns() - 1;

	int sampleInterval[2] = {1, 1};
	double var = CubicBezierSurfaceUtils::GetHeightVariance(grid, rowRange, colRange, sampleInterval);
	root->SetVariance(var);

	double maxVariance = var / resolution;
	// Make quadtree
	int currentDepth = 0;
	int maxDepth = 0;
	int nSamples = grid.size();
	bool resumeSubdivision = resolution <= 1 ? false : true;
	while (resumeSubdivision){
		resumeSubdivision = false;
		QuadTree::QuadTreeIterator iter((QuadTree *)root);
		currentDepth++;
		depthRange[0] = std::numeric_limits<int>::max();
		depthRange[1] = currentDepth;
		while ((*iter) != NULL){
			BezierQuadTree * bTree = (BezierQuadTree *)*iter;
			if (bTree->HaveChildren() == true || iter.GetDepth() >= currentDepth){
				++iter;
				continue;
			}

			rowRange = bTree->GetRowRange();
			colRange = bTree->GetColRange();

			var = CubicBezierSurfaceUtils::GetHeightVariance(grid, rowRange, colRange, sampleInterval);
			bTree->SetVariance(var);

			bool doMakeChildren = false;

			if (var > maxVariance){
				doMakeChildren = true;
				resumeSubdivision = true;
			}

			if (doMakeChildren && rowRange[1] - rowRange[0] + 1 >= minGridPoints && colRange[1] - colRange[0] + 1 >= minGridPoints){
				bTree->MakeChildren();

				maxDepth = currentDepth > maxDepth ? currentDepth : maxDepth;
				nSamples = rowRange[1] - rowRange[0] + 1 < nSamples ? rowRange[1] - rowRange[0] + 1 : nSamples;
				nSamples = colRange[1] - colRange[0] + 1 < nSamples ? colRange[1] - colRange[0] + 1 : nSamples;

				//// Make sure depth of neighbor is within -+1 difference
				QuadTree::QuadTreeIterator subIter((QuadTree *)root);
				if ((*subIter.GoTo(rowRange[0] - 1, colRange[0] + 1)) != NULL){
					BezierQuadTree * subBTree = (BezierQuadTree *)*subIter;
					if (iter.GetDepth() - subIter.GetDepth() > 0){
						subBTree->MakeChildren();
					}
				}
				subIter = root;
				if ((*subIter.GoTo(rowRange[1] + 1, colRange[0] + 1)) != NULL){
					BezierQuadTree * subBTree = (BezierQuadTree *)*subIter;
					if (iter.GetDepth() - subIter.GetDepth() > 0){
						subBTree->MakeChildren();
					}
				}
				subIter = root;
				if ((*subIter.GoTo(rowRange[0] + 1, colRange[0] - 1)) != NULL){
					BezierQuadTree * subBTree = (BezierQuadTree *)*subIter;
					if (iter.GetDepth() - subIter.GetDepth() > 0){
						subBTree->MakeChildren();
					}
				}
				subIter = root;
				if ((*subIter.GoTo(rowRange[0] + 1, colRange[1] + 1)) != NULL){
					BezierQuadTree * subBTree = (BezierQuadTree *)*subIter;
					if (iter.GetDepth() - subIter.GetDepth() > 0){
						subBTree->MakeChildren();
					}
				}
			}
			else{
				resumeSubdivision = false;
			}

			depthRange[0] = depthRange[0] > iter.GetDepth() ? iter.GetDepth() : depthRange[0];
			++iter;
		}
	}

	int nPatches = 0;
	QuadTree::QuadTreeIterator iter((QuadTree *)root);
	while ((*iter) != NULL){
		BezierQuadTree * bTree = (BezierQuadTree *)*iter;
		if (bTree->HaveChildren() == false){
			nPatches++;
		}
		++iter;
	}

	cout << "Made quad tree." << std::endl
		<< " Max depth:" << currentDepth << std::endl
		<< " Min. number of samples: " << nSamples << std::endl
		<< " Total # patches: " << nPatches << std::endl; 

	return root;
}

void
CubicBezierSurfaceUtils::fitPatchToGrid(BezierQuadTree::QuadTreeIterator & iter, Grid2D & grid, std::set< int > & gridIndexToSmooth, int shallowestDepth)
{
	BezierQuadTree * bTree = (BezierQuadTree *)*iter;

	BezierQuadTree * root = (BezierQuadTree *)(*iter)->GetRoot();
	set< int > depths;

	CubicBezierSurface * bezierSurface = bTree->GetSurface();
	Vector3D controlPoints[16];
	for (int i = 0; i < 16; i++){
		controlPoints[i] = bezierSurface->GetControlPoints()[i];
	}
	bool revert[16] = {
		false, false, false, false,
		false, false, false, false,
		false, false, false, false,
		false, false, false, false
	};

	int * rowRange = bTree->GetRowRange();
	int * colRange = bTree->GetColRange();
	int nRows = rowRange[1] - rowRange[0] + 1;
	int nCols = colRange[1] - colRange[0] + 1;
	int interval[2] = {(int)max(1.0, (double)nRows / 16.0), (int)max(1.0, (double)nCols / 16.0)};

	// Control points arrangement
	// 01           11 (row = 1)
	//   12 13 14 15
	//   08 09 10 11
	//   04 05 06 07
	//   00 01 02 03
	// 00           10 (row = 0)
	depths = getDepthsAroundGrid(grid, root, rowRange[0], colRange[0] + 2, iter.GetDepth(), -1);
	if (depths.size() < 2){
		FitBoundaryToGrid(*bezierSurface, grid, 0, 0, rowRange, colRange, interval[0]);
	}
	else{
		revert[0] = true;
		revert[3] = true;

		revert[4] = true;
		revert[5] = true;
		revert[6] = true;
		revert[7] = true;
	}

	depths = getDepthsAroundGrid(grid, root, rowRange[1], colRange[0] + 2, iter.GetDepth(), -1);
	if (depths.size() < 2){
		FitBoundaryToGrid(*bezierSurface, grid, 0, 1, rowRange, colRange, interval[0]);
	}
	else{
		revert[12] = true;
		revert[15] = true;

		revert[8] = true;
		revert[9] = true;
		revert[10] = true;
		revert[11] = true;
	}

	depths = getDepthsAroundGrid(grid, root, rowRange[0] + 2, colRange[0], iter.GetDepth(), -1);
	if (depths.size() < 2){
		FitBoundaryToGrid(*bezierSurface, grid, 1, 0, rowRange, colRange, interval[1]);
	}
	else{
		revert[0] = true;
		revert[12] = true;

		revert[1] = true;
		revert[5] = true;
		revert[9] = true;
		revert[13] = true;
	}

	depths = getDepthsAroundGrid(grid, root, rowRange[0] + 2, colRange[1], iter.GetDepth(), -1);
	if (depths.size() < 2){
		FitBoundaryToGrid(*bezierSurface, grid, 1, 1, rowRange, colRange, interval[1]);
	}
	else{
		revert[3] = true;
		revert[15] = true;

		revert[2] = true;
		revert[6] = true;
		revert[10] = true;
		revert[14] = true;
	}

	// corners
	depths = getDepthsAroundGrid(grid, root, rowRange[0], colRange[0], iter.GetDepth(), 0);
	if (depths.size() > 1){
		revert[0] = true;
		revert[1] = true;
		revert[4] = true;
		revert[5] = true;
	}
	depths = getDepthsAroundGrid(grid, root, rowRange[0], colRange[1], iter.GetDepth(), 0);
	if (depths.size() > 1){
		revert[2] = true;
		revert[3] = true;
		revert[6] = true;
		revert[7] = true;
	}
	depths = getDepthsAroundGrid(grid, root, rowRange[1], colRange[0], iter.GetDepth(), 0);
	if (depths.size() > 1){
		revert[8] = true;
		revert[9] = true;
		revert[12] = true;
		revert[13] = true;
	}
	depths = getDepthsAroundGrid(grid, root, rowRange[1], colRange[1], iter.GetDepth(), 0);
	if (depths.size() > 1){
		revert[10] = true;
		revert[11] = true;
		revert[14] = true;
		revert[15] = true;
	}

	AdjustInternalControlPointsToFitGrid(*bezierSurface, grid, rowRange, colRange, interval);

	if (iter.GetDepth() > shallowestDepth){
		for (int i = 0; i < 16; i++){
			if (revert[i]){
				bezierSurface->GetControlPoints()[i] = controlPoints[i];
			}

			if (i == 0 && !revert[i]){
				gridIndexToSmooth.insert(grid.getIndex(colRange[0], rowRange[0]));
			}
			if (i == 3 && !revert[i]){
				gridIndexToSmooth.insert(grid.getIndex(colRange[1], rowRange[0]));
			}
			if (i == 12 && !revert[i]){
				gridIndexToSmooth.insert(grid.getIndex(colRange[0], rowRange[1]));
			}
			if (i == 15 && !revert[i]){
				gridIndexToSmooth.insert(grid.getIndex(colRange[1], rowRange[1]));
			}
		}
	}
}


void
CubicBezierSurfaceUtils::doGlobalSmoothing(BezierQuadTree * root, Grid2D &grid, set< int > & gridIndexToSmooth, int depthLimit)
{
	// Smooth to achieve G1 continuity
	std::map<int, vector < vector< Vector3D * > > > nodeToControlPoints; // grid node to neighbor control points reference
	BezierQuadTree::QuadTreeIterator iter(root);

	int rootDepth = iter.GetDepth();
	int upperLimit = depthLimit == -1 ? -1 : rootDepth + depthLimit;
	int lowerLimit = depthLimit == -1 ? -1 : rootDepth;

	++iter;
	while(iter != NULL){
		int currentDepth = iter.GetDepth();
		BezierQuadTree * bTree = (BezierQuadTree *)*iter;

		if ( (lowerLimit != -1 && currentDepth < lowerLimit) || (upperLimit != -1 && upperLimit < currentDepth) ){
			++iter;
			continue;
		}

		if (bTree->HaveChildren() && currentDepth < upperLimit ){
			++iter;
			continue;
		}

		if (bTree->HaveSurface() == false){
			++iter;
			continue;
		}

		int * rowRange = bTree->GetRowRange();
		int * colRange = bTree->GetColRange();

		int index;
		vector< vector< Vector3D * > > points;
		CubicBezierSurface * bezierSurface = bTree->GetSurface();
		Vector3D * controlPoints = bezierSurface->GetControlPoints();

		// Control points arrangement
		// 01           11 (row = 1)
		//   12 13 14 15
		//   08 09 10 11
		//   04 05 06 07
		//   00 01 02 03
		// 00           10 (row = 0)

		// Sort the points
		// 6 7 8
		// 3 4 5
		// 0 1 2

		// Corner col=0, row=0
		index = grid.getIndex(colRange[0], rowRange[0]);
		set< int > depths = getDepthsAroundGrid(grid, root, rowRange[0], colRange[0], upperLimit, lowerLimit);
		if ( (gridIndexToSmooth.size() == 0 || gridIndexToSmooth.find( index ) != gridIndexToSmooth.end()) && depths.size() == 1){
			points = vector < vector< Vector3D *> >();
			points.resize(9);
			if (nodeToControlPoints.find(index) != nodeToControlPoints.end()){
				points = nodeToControlPoints[index];
			}
			points[4].push_back(&(controlPoints[0]));
			points[5].push_back(&(controlPoints[1]));
			points[7].push_back(&(controlPoints[4]));
			points[8].push_back(&(controlPoints[5]));
			nodeToControlPoints[index] = points;
		}

		// Corner col=1, row=0
		index = grid.getIndex(colRange[1], rowRange[0]);
		depths = getDepthsAroundGrid(grid, root, rowRange[0], colRange[1], upperLimit, lowerLimit);
		if ( (gridIndexToSmooth.size() == 0 || gridIndexToSmooth.find( index ) != gridIndexToSmooth.end()) && depths.size() == 1){
			points = vector < vector< Vector3D *> >();
			points.resize(9);
			if (nodeToControlPoints.find(index) != nodeToControlPoints.end()){
				points = nodeToControlPoints[index];
			}
			points[3].push_back(&(controlPoints[2]));
			points[4].push_back(&(controlPoints[3]));
			points[6].push_back(&(controlPoints[6]));
			points[7].push_back(&(controlPoints[7]));
			nodeToControlPoints[index] = points;
		}

		// Corner col=0, row=1
		index = grid.getIndex(colRange[0], rowRange[1]);
		depths = getDepthsAroundGrid(grid, root, rowRange[1], colRange[0], upperLimit, lowerLimit);
		if ( (gridIndexToSmooth.size() == 0 || gridIndexToSmooth.find( index ) != gridIndexToSmooth.end()) && depths.size() == 1){
			points = vector < vector< Vector3D *> >();
			points.resize(9);
			if (nodeToControlPoints.find(index) != nodeToControlPoints.end()){
				points = nodeToControlPoints[index];
			}
			points[1].push_back(&(controlPoints[8]));
			points[2].push_back(&(controlPoints[9]));
			points[4].push_back(&(controlPoints[12]));
			points[5].push_back(&(controlPoints[13]));
			nodeToControlPoints[index] = points;
		}

		// Corner col=1, row=1
		index = grid.getIndex(colRange[1], rowRange[1]);
		depths = getDepthsAroundGrid(grid, root, rowRange[1], colRange[1], upperLimit, lowerLimit);
		if ( (gridIndexToSmooth.size() == 0 || gridIndexToSmooth.find( index ) != gridIndexToSmooth.end()) && depths.size() == 1){
			points = vector < vector< Vector3D *> >();
			points.resize(9);
			if (nodeToControlPoints.find(index) != nodeToControlPoints.end()){
				points = nodeToControlPoints[index];
			}
			points[0].push_back(&(controlPoints[10]));
			points[1].push_back(&(controlPoints[11]));
			points[3].push_back(&(controlPoints[14]));
			points[4].push_back(&(controlPoints[15]));
			nodeToControlPoints[index] = points;
		}

		++iter;
	}

	smoothPoints(nodeToControlPoints);
}

void
CubicBezierSurfaceUtils::alignAllPatchesToGrid(BezierQuadTree * root, int depth)
{
	BezierQuadTree::QuadTreeIterator iter(root);
	CubicBezierSurface * bezierSurface = root->GetSurface();
	PSpace pSpace = bezierSurface->GetPSpace();

	++iter;
	while(iter != NULL){
		BezierQuadTree * bTree = (BezierQuadTree *)*iter;
		if (!bTree->HaveSurface()){
			++iter;
			continue;
		}
		if (depth != -1 && iter.GetDepth() != depth){
			continue;
		}
		alignPatchToGrid(bTree);

		++iter;
	}
}

void
CubicBezierSurfaceUtils::alignPatchToGrid(BezierQuadTree * bTree)
{
	if (!bTree->HaveSurface()){
		return;
	}

	int * rowRange = bTree->GetRowRange();
	int * colRange = bTree->GetColRange();
	CubicBezierSurface * bezierSurface = bTree->GetSurface();
	PSpace pSpace = bezierSurface->GetPSpace();
	double xRange[2];
	double yRange[2];
	for (int i = 0; i < 2; i++){
		xRange[i] = colRange[i] * pSpace.GetAxisLength(0);
		yRange[i] = rowRange[i] * pSpace.GetAxisLength(1);
	}
	for (int i = 0; i < 4; i++){
		bezierSurface->GetControlPoints()[4 * i][0] = xRange[0];
		bezierSurface->GetControlPoints()[4 * i + 3][0] = xRange[1];
		bezierSurface->GetControlPoints()[     i][1] = yRange[0];
		bezierSurface->GetControlPoints()[12 + i][1] = yRange[1];
	}
}

void CubicBezierSurfaceUtils::smoothPoints(std::map<int,  std::vector< std::vector < Vector3D * > > > & nodeToControlPoints)
{
	bool doSmoothest = true;

	for (map<int,  vector< vector < Vector3D * > > >::iterator nodeToControlPointsIter(nodeToControlPoints.begin()); nodeToControlPointsIter != nodeToControlPoints.end(); nodeToControlPointsIter++){
		vector< vector< Vector3D * > > sortedPoints = nodeToControlPointsIter->second;
		vector< Vector3D * > points;
		for (int i = 0; i < (int)sortedPoints.size(); i++){
			for (int j = 0; j < (int)sortedPoints[i].size(); j++){
				points.push_back( sortedPoints[i][j] );
			}
		}
		if (points.size() <= 4){
			continue;
		}
		Vector3D gridPoint = *(sortedPoints[4][0]);

		Matrix4X4 m2;
		Vector3D z;
		for (int i = 0; i < (int)points.size(); i++){
			Vector3D coefficient((*points[i])[0], (*points[i])[1], 1);
			for (int j = 0; j < 3; j++){
				m2[j][0] += coefficient[0] * coefficient[j];
				m2[j][1] += coefficient[1] * coefficient[j];
				m2[j][2] += coefficient[2] * coefficient[j];
				z[j] -= (*points[i])[2] * coefficient[j];
			}
		}
		m2[3][3] = 1.0;
		Matrix4X4 invM2 = invertMatrix4X4( m2 );

		Vector3D norm(0.0, 0.0, 1);
		double d = 0.0;
		for (int i = 0; i < 3; i++){
			norm[0] += invM2[0][i] * z[i];
			norm[1] += invM2[1][i] * z[i];
			d += invM2[2][i] * z[i];
		}

		int index[5] = {4, 1, 3, 7, 5};
		for (int i = 0; i < 5; i++){
			if (sortedPoints[index[i]].size() == 0){
				continue;
			}
			double x0 = (*(sortedPoints[index[i]][0]))[0];
			double y0 = (*(sortedPoints[index[i]][0]))[1];
			(*sortedPoints[index[i]][0])[2] = -(norm[0] * x0 + norm[1] * y0 + d) / norm[2];
			for (int j = 1; j < (int)sortedPoints[index[i]].size(); j++){
				(*(sortedPoints[index[i]][j])) = (*(sortedPoints[index[i]][0]));
			}
		}
		// 6 7 8
		// 3 4 5
		// 0 1 2
		// Among sortedPoints[0, 2, 6, 8], find the closest height to sortedPoints[4]
		index[0] = 0;
		index[1] = 2;
		index[2] = 6;
		index[3] = 8;
		map<double, int> order;
		for (int i = 0; i < 4; i++){
			if (sortedPoints[index[i]].size() > 0){
				order[ abs((*sortedPoints[index[i]][0])[2] - (*sortedPoints[4][0])[2]) ] = index[i];
			}
		}
		Vector3D u;
		Vector3D v;
		switch((*order.begin()).second){
			case 0:
				u = (*sortedPoints[3][0]) - (*sortedPoints[4][0]);
				v = (*sortedPoints[1][0]) - (*sortedPoints[4][0]);
				for (int i = 0; i < (int) sortedPoints[0].size(); i++){
					(*sortedPoints[0][i]) = (*sortedPoints[4][0]) + u + v;
				}
				for (int i = 0; i < (int) sortedPoints[6].size(); i++){
					(*sortedPoints[6][i]) = (*sortedPoints[4][0]) + u - v;
				}
				for (int i = 0; i < (int) sortedPoints[7].size(); i++){
					(*sortedPoints[7][i]) = (*sortedPoints[4][0]) - v;
				}
				for (int i = 0; i < (int) sortedPoints[8].size(); i++){
					(*sortedPoints[8][i]) = (*sortedPoints[4][0]) - u - v;
				}
				for (int i = 0; i < (int) sortedPoints[5].size(); i++){
					(*sortedPoints[5][i]) = (*sortedPoints[4][0]) - u;
				}
				for (int i = 0; i < (int) sortedPoints[2].size(); i++){
					(*sortedPoints[2][i]) = (*sortedPoints[4][0]) - u + v;
				}
				break;
			case 6:
				u = (*sortedPoints[3][0]) - (*sortedPoints[4][0]);
				v = (*sortedPoints[7][0]) - (*sortedPoints[4][0]);
				for (int i = 0; i < (int) sortedPoints[6].size(); i++){
					(*sortedPoints[6][i]) = (*sortedPoints[4][0]) + u + v;
				}
				for (int i = 0; i < (int) sortedPoints[8].size(); i++){
					(*sortedPoints[8][i]) = (*sortedPoints[4][0]) - u + v;
				}
				for (int i = 0; i < (int) sortedPoints[5].size(); i++){
					(*sortedPoints[5][i]) = (*sortedPoints[4][0]) - u;
				}
				for (int i = 0; i < (int) sortedPoints[2].size(); i++){
					(*sortedPoints[2][i]) = (*sortedPoints[4][0]) - u - v;
				}
				for (int i = 0; i < (int) sortedPoints[1].size(); i++){
					(*sortedPoints[1][i]) = (*sortedPoints[4][0]) - v;
				}
				for (int i = 0; i < (int) sortedPoints[0].size(); i++){
					(*sortedPoints[0][i]) = (*sortedPoints[4][0]) + u - v;
				}
				break;
			case 8:
				u = (*sortedPoints[5][0]) - (*sortedPoints[4][0]);
				v = (*sortedPoints[7][0]) - (*sortedPoints[4][0]);
				for (int i = 0; i < (int) sortedPoints[8].size(); i++){
					(*sortedPoints[8][i]) = (*sortedPoints[4][0]) + u + v;
				}
				for (int i = 0; i < (int) sortedPoints[2].size(); i++){
					(*sortedPoints[2][i]) = (*sortedPoints[4][0]) + u - v;
				}
				for (int i = 0; i < (int) sortedPoints[1].size(); i++){
					(*sortedPoints[1][i]) = (*sortedPoints[4][0]) - v;
				}
				for (int i = 0; i < (int) sortedPoints[0].size(); i++){
					(*sortedPoints[0][i]) = (*sortedPoints[4][0]) - u - v;
				}
				for (int i = 0; i < (int) sortedPoints[3].size(); i++){
					(*sortedPoints[3][i]) = (*sortedPoints[4][0]) - u;
				}
				for (int i = 0; i < (int) sortedPoints[6].size(); i++){
					(*sortedPoints[6][i]) = (*sortedPoints[4][0]) - u + v;
				}
				break;
		// 6 7 8
		// 3 4 5
		// 0 1 2
			case 2:
				u = (*sortedPoints[5][0]) - (*sortedPoints[4][0]);
				v = (*sortedPoints[1][0]) - (*sortedPoints[4][0]);
				for (int i = 0; i < (int) sortedPoints[2].size(); i++){
					(*sortedPoints[2][i]) = (*sortedPoints[4][0]) + u + v;
				}
				for (int i = 0; i < (int) sortedPoints[0].size(); i++){
					(*sortedPoints[0][i]) = (*sortedPoints[4][0]) - u + v;
				}
				for (int i = 0; i < (int) sortedPoints[3].size(); i++){
					(*sortedPoints[3][i]) = (*sortedPoints[4][0]) - u;
				}
				for (int i = 0; i < (int) sortedPoints[6].size(); i++){
					(*sortedPoints[6][i]) = (*sortedPoints[4][0]) - u - v;
				}
				for (int i = 0; i < (int) sortedPoints[7].size(); i++){
					(*sortedPoints[7][i]) = (*sortedPoints[4][0]) - v;
				}
				for (int i = 0; i < (int) sortedPoints[8].size(); i++){
					(*sortedPoints[8][i]) = (*sortedPoints[4][0]) + u - v;
				}
				break;
		}
	}
}


void
CubicBezierSurfaceUtils::makeChildSubdivisionSurface(BezierQuadTree * bTree)
{
	bool fineTuneParam = true;
	bool adjustBoundary = true;
	if (!bTree->HaveChildren()){
		return;
	}

	CubicBezierSurface * bezierSurface = bTree->GetSurface();
	int * parentRowRange;
	int * parentColRange;
	parentRowRange = bTree->GetRowRange();
	parentColRange = bTree->GetColRange();
	Vector3D * parentControlPoints = bezierSurface->GetControlPoints();

	BezierQuadTree * childTree;
	double uRange[2];
	double vRange[2];
	int * rowRange;
	int * colRange;
	Vector3D * controlPoints;

	int index[2];

	// ^
	// |
	// r
	// o
	// w
	// * 4 *
	// 2 0 3
	// * 1 *
	// col->
	double param[5][2];


	if (fineTuneParam){
		childTree = (BezierQuadTree * )bTree->GetChildren()[0];
		rowRange = childTree->GetRowRange();
		colRange = childTree->GetColRange();
		bTree->GridToParameters((double)rowRange[1], (double)colRange[1], index, param[0], 0);
		bTree->GridToParameters((double)parentRowRange[0], (double)colRange[1], index, param[1], 0); // (row, col) = (0.0, 0.5)
		bTree->GridToParameters((double)rowRange[1], (double)parentColRange[0], index, param[2], 0); // (row, col) = (0.5, 0.0)
		bTree->GridToParameters((double)rowRange[1], (double)parentColRange[1], index, param[3], 0); // (row, col) = (0.5, 1.0)
		bTree->GridToParameters((double)parentRowRange[1], (double)colRange[1], index, param[4], 0); // (row, col) = (1.0, 0.5)
	}

	// 1 3
	// 0 2
	// Child 0
	uRange[0] = 0.0;
	uRange[1] = param[0][0];
	vRange[0] = 0.0;
	vRange[1] = param[0][1];
	childTree = (BezierQuadTree * )bTree->GetChildren()[0];
	childTree->SetSurface( new CubicBezierSurface( CubicBezierSurfaceUtils::Divide(*bezierSurface, uRange, vRange) ) );
	if (adjustBoundary){
		controlPoints = childTree->GetSurface()->GetControlPoints();
		CubicBezierCurve3D bezierCurve( bTree->GetSurface()->GetBoundaryCurve(0) );
		bezierCurve = CubicBezierCurve3DUtils::Divide(bezierCurve, 0.0, param[1][0]);
		for (int i = 0; i < 4; i++){
			controlPoints[i] = bezierCurve.GetControlPoints()[i];
		}
		bezierCurve = bTree->GetSurface()->GetBoundaryCurve(2);
		bezierCurve = CubicBezierCurve3DUtils::Divide(bezierCurve, 0.0, param[2][1]);
		for (int i = 0; i < 4; i++){
			controlPoints[4 * i] = bezierCurve.GetControlPoints()[i];
		}
	}

	// Child 1
	uRange[0] = param[0][0];
	uRange[1] = 1.0;
	vRange[0] = 0.0;
	vRange[1] = param[0][1];
	childTree = (BezierQuadTree * )bTree->GetChildren()[1];
	childTree->SetSurface( new CubicBezierSurface( CubicBezierSurfaceUtils::Divide(*bezierSurface, uRange, vRange) ) );
	if (adjustBoundary){
		controlPoints = childTree->GetSurface()->GetControlPoints();
		CubicBezierCurve3D bezierCurve( bTree->GetSurface()->GetBoundaryCurve(0) );
		bezierCurve = CubicBezierCurve3DUtils::Divide(bezierCurve, param[1][0], 1.0);
		for (int i = 0; i < 4; i++){
			controlPoints[i] = bezierCurve.GetControlPoints()[i];
		}
		bezierCurve = bTree->GetSurface()->GetBoundaryCurve(3);
		bezierCurve = CubicBezierCurve3DUtils::Divide(bezierCurve, 0.0, param[3][1]);
		for (int i = 0; i < 4; i++){
			controlPoints[4 * i + 3] = bezierCurve.GetControlPoints()[i];
		}
	}

	// Child 2
	uRange[0] = 0.0;
	uRange[1] = param[0][0];
	vRange[0] = param[0][1];
	vRange[1] = 1.0;
	childTree = (BezierQuadTree * )bTree->GetChildren()[2];
	childTree->SetSurface( new CubicBezierSurface( CubicBezierSurfaceUtils::Divide(*bezierSurface, uRange, vRange) ) );
	if (adjustBoundary){
		controlPoints = childTree->GetSurface()->GetControlPoints();
		CubicBezierCurve3D bezierCurve( bTree->GetSurface()->GetBoundaryCurve(1) );
		bezierCurve = CubicBezierCurve3DUtils::Divide(bezierCurve, 0.0, param[4][0]);
		for (int i = 0; i < 4; i++){
			controlPoints[12 + i] = bezierCurve.GetControlPoints()[i];
		}
		bezierCurve = bTree->GetSurface()->GetBoundaryCurve(2);
		bezierCurve = CubicBezierCurve3DUtils::Divide(bezierCurve, param[2][1], 1.0);
		for (int i = 0; i < 4; i++){
			controlPoints[4 * i] = bezierCurve.GetControlPoints()[i];
		}
	}

	// Child 3
	uRange[0] = param[0][0];
	uRange[1] = 1.0;
	vRange[0] = param[0][1];
	vRange[1] = 1.0;
	childTree = (BezierQuadTree * )bTree->GetChildren()[3];
	childTree->SetSurface( new CubicBezierSurface( CubicBezierSurfaceUtils::Divide(*bezierSurface, uRange, vRange) ) );
	if (adjustBoundary){
		controlPoints = childTree->GetSurface()->GetControlPoints();
		CubicBezierCurve3D bezierCurve( bTree->GetSurface()->GetBoundaryCurve(1) );
		bezierCurve = CubicBezierCurve3DUtils::Divide(bezierCurve, param[4][0], 1.0);
		for (int i = 0; i < 4; i++){
			controlPoints[12 + i] = bezierCurve.GetControlPoints()[i];
		}
		bezierCurve = bTree->GetSurface()->GetBoundaryCurve(3);
		bezierCurve = CubicBezierCurve3DUtils::Divide(bezierCurve, param[3][1], 1.0);
		for (int i = 0; i < 4; i++){
			controlPoints[4 * i + 3] = bezierCurve.GetControlPoints()[i];
		}
	}
}

set< int >
CubicBezierSurfaceUtils::getDepthsAroundGrid(Grid2D &grid, BezierQuadTree * root, int row, int col, int upperLimit, int lowerLimit)
{
	set< int > depths;
	upperLimit = upperLimit == -1 ? std::numeric_limits<int>::max() : upperLimit;
	lowerLimit = lowerLimit == -1 ? -std::numeric_limits<int>::max() : lowerLimit;

	for (int i = -1; i <= 1; i++){
		for (int j = -1; j <= 1; j++){
			BezierQuadTree::QuadTreeIterator iter(root);
			iter.GoTo(row + i, col + j);
			if (iter != NULL){
				int depth = iter.GetDepth();
				if (lowerLimit <= depth && depth <= upperLimit){
					depths.insert( iter.GetDepth() );
				}
			}
		}
	}

	return depths;
}

static std::set< CubicBezierSurface * >getBezierSurfacesAroundGrid(Grid2D &grid, BezierQuadTree * root, int row, int col)
{
	set< CubicBezierSurface * > surfaces;

	for (int i = -1; i <= 1; i++){
		for (int j = -1; j <= 1; j++){
			BezierQuadTree::QuadTreeIterator iter(root);
			iter.GoTo(row + i, col + j);
			if (iter != NULL){
				BezierQuadTree * bTree = (BezierQuadTree *)*iter;
				surfaces.insert( bTree->GetSurface() );
			}
		}
	}

	return surfaces;
}

CubicBezierSurface
CubicBezierSurfaceUtils::FitSurfaceToGrid(Grid2D & grid, int rowRange[2], int colRange[2], int sampleInterval[2])
{
	if (sampleInterval[0] == 0 || sampleInterval[1] == 0){
		return NULL;
	}
	if (rowRange[0] >= rowRange[1] || rowRange[0] < 0 || (int)grid.numRows() <= rowRange[1]){
		return NULL;
	}
	if (colRange[0] >= colRange[1] || colRange[0] < 0 || (int)grid.numColumns() <= colRange[1]){
		return NULL;
	}

	Vector3D cp[16];
	CubicBezierSurface bezierSurface(cp);
	bezierSurface.SetPSpace( grid.GetPSpace() );
	Vector3D * controlPoints = bezierSurface.GetControlPoints();

	for (int i = 0; i < 16; i++){
		controlPoints[i] = Vector3D(0.0, 0.0, 0.0);
	}

	// First, get the 4 boundary Bezier curves. Storing u0, u1, v0, v1 order
	// u
	for (int row = 0; row < 2; row++){
		FitBoundaryToGrid(bezierSurface, grid, 0, row, rowRange, colRange, sampleInterval[0]);
	}
	// v
	for (int col = 0; col < 2; col++){
		FitBoundaryToGrid(bezierSurface, grid, 1, col, rowRange, colRange, sampleInterval[1]);
	}

	AdjustInternalControlPointsToFitGrid(bezierSurface, grid, rowRange, colRange, sampleInterval);

	return bezierSurface;
}

void
CubicBezierSurfaceUtils::FitBoundaryToGrid(CubicBezierSurface & bezier, Grid2D & grid, int dir, int n, int rowRange[2], int colRange[2], int sampleInterval)
{
	if (sampleInterval == 0){
		return;
	}
	if (rowRange[0] >= rowRange[1] || rowRange[0] < 0 || (int)grid.numRows() <= rowRange[1]){
		return;
	}
	if (colRange[0] >= colRange[1] || colRange[0] < 0 || (int)grid.numColumns() <= colRange[1]){
		return;
	}

	Vector3D * controlPoints = bezier.GetControlPoints();

	vector<Vector3D> polyline;
	if (dir == 0){
		// u direction corresponds to column direction
		int row = n;
		ExtractGridPolyline(grid, 0, rowRange[row], colRange, sampleInterval, polyline);
		CubicBezierCurve3D bezier = CubicBezierCurve3DUtils::FitCurveToPolyline(polyline);
		for (int i = 0; i < 4; i++){
			controlPoints[12 * row + i] = bezier.GetControlPoints()[i];
		}
	}
	else{
		// v direction corresponds to row direction
		int col = n;
		ExtractGridPolyline(grid, 1, colRange[col], rowRange, sampleInterval, polyline);
		CubicBezierCurve3D bezier = CubicBezierCurve3DUtils::FitCurveToPolyline(polyline);
		for (int i = 0; i < 4; i++){
			controlPoints[4 * i + 3 * col ] = bezier.GetControlPoints()[i];
		}
	}
}

void
CubicBezierSurfaceUtils::AdjustInternalControlPointsToFitGrid(CubicBezierSurface & bezier, Grid2D & grid, int rowRange[2], int colRange[2], int sampleInterval[2])
{
	Vector3D * controlPoints = bezier.GetControlPoints();

	if (sampleInterval[0] == 0 || sampleInterval[1] == 0){
		return;
	}
	if (rowRange[0] >= rowRange[1] || rowRange[0] < 0 || (int)grid.numRows() <= rowRange[1]){
		return;
	}
	if (colRange[0] >= colRange[1] || colRange[0] < 0 || (int)grid.numColumns() <= colRange[1]){
		return;
	}

	double rangeX[2] = {MAX_FLOAT, -MAX_FLOAT};
	double rangeY[2] = {MAX_FLOAT, -MAX_FLOAT};
	double rangeZ[2] = {MAX_FLOAT, -MAX_FLOAT};
	for (int i = 0; i < 16; i++){
		if ( (i / 4 != 1 && i / 4 != 2) || (i % 4 != 1 && i % 4 != 2) ){
			rangeX[0] =  rangeX[0] > controlPoints[i][0] ? controlPoints[i][0] : rangeX[0];
			rangeX[1] =  rangeX[1] < controlPoints[i][0] ? controlPoints[i][0] : rangeX[1];
			rangeY[0] =  rangeY[0] > controlPoints[i][1] ? controlPoints[i][1] : rangeY[0];
			rangeY[1] =  rangeY[1] < controlPoints[i][1] ? controlPoints[i][1] : rangeY[1];
			rangeZ[0] =  rangeZ[0] > controlPoints[i][2] ? controlPoints[i][2] : rangeZ[0];
			rangeZ[1] =  rangeZ[1] < controlPoints[i][2] ? controlPoints[i][2] : rangeZ[1];
		}
	}

	for (int i = 1; i < 3; i++){
		for (int j = 1; j < 3; j++){
			controlPoints[i * 4 + j] = Vector3D(0.0, 0.0, 0.0);
		}
	}

	int nSamplesRow = (int)ceil( (double)(rowRange[1] - rowRange[0] + 1) / (double)sampleInterval[0]);
	int nSamplesCol = (int)ceil( (double)(colRange[1] - colRange[0] + 1) / (double)sampleInterval[1]);

	// Compute the parameters for the sample point
	vector<Vector3D> polyline;
	vector< double > u;
	u.resize( (nSamplesRow - 2) * (nSamplesCol - 2) );
	vector< double > v;
	v.resize( u.size() );
	vector< Vector3D > samples;
	samples.resize( u.size() );

	for (int row = 0; row < (nSamplesRow - 2); row++){
		ExtractGridPolyline(grid, 0, rowRange[0] + 1 + row * sampleInterval[0], colRange, sampleInterval[1], polyline);
		for (int col = 0; col < (nSamplesCol - 2); col++){
			samples[row * (nSamplesCol - 2) + col] = grid.getPoint( grid.getIndex(colRange[0] + 1 + col * sampleInterval[1], rowRange[0] + 1 + row * sampleInterval[0] ) );
			u[row * (nSamplesCol - 2) + col] = CubicBezierCurve3DUtils::PolylinePointToParam(polyline, col + 1, true);
		}
	}
	for (int col = 0; col < (nSamplesCol - 2); col++){
		ExtractGridPolyline(grid, 1, colRange[0] + 1 + col * sampleInterval[1], rowRange, sampleInterval[0], polyline);
		for (int row = 0; row < (nSamplesRow - 2); row++){
			v[row * (nSamplesCol - 2) + col] = CubicBezierCurve3DUtils::PolylinePointToParam(polyline, row + 1, true);
		}
	}

	// Compute the known q(k)
	vector< Vector3D > q;
	q.resize( u.size() );
	// Build 4x4 matrix of Bernstein bm = [bubv1, bubv2, bubv3, bubv4]TT x [bubv1, bubv2, bubv3, bubv4]T
	Matrix4X4 bm;
	// r = sum of [bubv0, bubv1, bubv2, bubv3]TT * [q0, q1, q2, q3]T
	Vector3D r[4];
	CubicBezierSurface nilInnerBezier(controlPoints); // Bicubic Bezier function with inner control points being nil
	for (int i = 0; i < (int)q.size(); i++){
		q[i] = samples[i] - nilInnerBezier.GetPoint(u[i], v[i]);

		// bubv: Bernstein combinations coeffients for inner control points
		double bubv[4];
		double bu[2] = {
			CubicBezierCurve3D::Bernstein(u[i], 1),
			CubicBezierCurve3D::Bernstein(u[i], 2),
		};
		double bv[2] = {
			CubicBezierCurve3D::Bernstein(v[i], 1),
			CubicBezierCurve3D::Bernstein(v[i], 2),
		};
		bubv[0] = bu[0] * bv[0];
		bubv[1] = bu[0] * bv[1];
		bubv[2] = bu[1] * bv[0];
		bubv[3] = bu[1] * bv[1];

		for (int j = 0; j < 4; j++){
			bm[j][0] += bubv[0] * bubv[j];
			bm[j][1] += bubv[1] * bubv[j];
			bm[j][2] += bubv[2] * bubv[j];
			bm[j][3] += bubv[3] * bubv[j];

			r[j] += bubv[j] * q[i];
		}

	}

	Matrix4X4 invBm = invertMatrix4X4( bm );

	Vector3D p[4];
	for (int i = 0; i < 4; i++){
		p[0] += invBm[0][i] * r[i];
		p[1] += invBm[1][i] * r[i];
		p[2] += invBm[2][i] * r[i];
		p[3] += invBm[3][i] * r[i];
	}

	for (int i = 0; i < 4; i++){
		p[i][0] = p[i][0] > rangeX[0] ? p[i][0] :  rangeX[0];
		p[i][0] = p[i][0] < rangeX[1] ? p[i][0] :  rangeX[1];
		p[i][1] = p[i][1] > rangeY[0] ? p[i][1] :  rangeY[0];
		p[i][1] = p[i][1] < rangeY[1] ? p[i][1] :  rangeY[1];
		//p[i][2] = p[i][2] > rangeZ[0] - cushion * (rangeZ[1] - rangeZ[0]) ? p[i][2] : rangeZ[0] - cushion * (rangeZ[1] - rangeZ[0]);
		//p[i][2] = p[i][2] < rangeZ[1] + cushion * (rangeZ[1] - rangeZ[0]) ? p[i][2] : rangeZ[1] + cushion * (rangeZ[1] - rangeZ[0]);
	}

	controlPoints[ 4 * 1 + 1 ] = p[0];
	controlPoints[ 4 * 1 + 2 ] = p[1];
	controlPoints[ 4 * 2 + 1 ] = p[2];
	controlPoints[ 4 * 2 + 2 ] = p[3];
}

Matrix4X4
CubicBezierSurfaceUtils::invertMatrix4X4(Matrix4X4 & m)
{
	//		Matrix4X4 invm = Inv(m);
	double det = Det( m );
	if (det == 0.0){
		throw Exception( "CubicBezierSurfaceUtils::invertMatrix4x4 determinant is zero." );
	}
	Matrix4X4 a(
		Minor(m, 1, 2, 3, 1, 2, 3),
		-Minor(m, 1, 2, 3, 0, 2, 3),
		Minor(m, 1, 2, 3, 0, 1, 3),
		Minor(m, 1, 2, 3, 0, 1, 2),
		-Minor(m, 0, 2, 3, 1, 2, 3),
		Minor(m, 0, 2, 3, 0, 2, 3),
		-Minor(m, 0, 2, 3, 0, 1, 3),
		Minor(m, 0, 2, 3, 0, 1, 2),
		Minor(m, 0, 1, 3, 1, 2, 3),
		-Minor(m, 0, 1, 3, 0, 2, 3),
		Minor(m, 0, 1, 3, 0, 1, 3),
		-Minor(m, 0, 1, 3, 0, 1, 2),
		-Minor(m, 0, 1, 2, 1, 2, 3),
		Minor(m, 0, 1, 2, 0, 2, 3),
		-Minor(m, 0, 1, 2, 0, 1, 3),
		Minor(m, 0, 1, 2, 0, 1, 2));
	Matrix4X4 invm = a / det;
	invm[3][0] = -invm[3][0];

	return invm;
}

void
CubicBezierSurfaceUtils::ExtractGridPolyline(Grid2D & grid, int dir, int n, int range[2], int interval, vector<Vector3D> & polyline)
{
	if (dir == 0){
		if (range[0] >= range[1] || range[0] < 0 || (int)grid.numColumns() <= range[1]) {
			return;
		}

		int row = n;
		int nSamplesCol = (int)ceil((double)(range[1] - range[0] + 1) / (double)interval);
		int col;
		polyline.resize( nSamplesCol );
		for (int i = 0; i < nSamplesCol - 1; i++){
			col = range[0] + interval * i;
			polyline[i] = grid.getPoint(grid.getIndex(col, row));
		}
		col = range[1];
		polyline[nSamplesCol - 1] = grid.getPoint(grid.getIndex(col, row));
	}
	else{
		if (range[0] >= range[1] || range[0] < 0 || (int)grid.numRows() <= range[1]) {
			return;
		}

		int col = n;
		int nSamplesRow = (int)ceil( (double)(range[1] - range[0] + 1) / (double)interval);
		int row;
		polyline.resize( nSamplesRow );
		for (int i = 0; i < nSamplesRow - 1; i++){
			row = range[0] + interval * i;
			polyline[i] = grid.getPoint(grid.getIndex(col, row));
		}
		row = range[1];
		polyline[nSamplesRow - 1] = grid.getPoint(grid.getIndex(col, row));
	}
}

double
CubicBezierSurfaceUtils::GetHeightVariance(Grid2D & grid, int rowRange[2], int colRange[2], int sampleInterval[2])
{
	if (sampleInterval[0] == 0 || sampleInterval[1] == 0){
		return NULL;
	}
	if (rowRange[0] >= rowRange[1] || rowRange[0] < 0 || (int)grid.numRows() <= rowRange[1]){
		return NULL;
	}
	if (colRange[0] >= colRange[1] || colRange[0] < 0 || (int)grid.numColumns() <= colRange[1]){
		return NULL;
	}

	int nSamplesRow = (int)ceil( (double)(rowRange[1] - rowRange[0] + 1) / (double)sampleInterval[0]);
	int nSamplesCol = (int)ceil( (double)(colRange[1] - colRange[0] + 1) / (double)sampleInterval[1]);
	int nSamples = nSamplesRow * nSamplesCol;

	double sumHeightDiff = 0;
	double maxHeight = - MAX_FLOAT;
	double minHeight = MAX_FLOAT;
	for (int row = rowRange[0] + 1; row < rowRange[1] - 1; row += sampleInterval[0]){
		for (int col = colRange[0] + 1; col < colRange[1] - 1; col+= sampleInterval[1]){
			double currentHeight = grid.getGridVal(grid.getIndex(col, row));

			double heightDiff = 0;
			heightDiff += abs(grid.getGridVal(grid.getIndex(col - 1, row - 1)) - currentHeight);
			heightDiff += abs(grid.getGridVal(grid.getIndex(col, row - 1)) - currentHeight);
			heightDiff += abs(grid.getGridVal(grid.getIndex(col + 1, row - 1)) - currentHeight);

			heightDiff += abs(grid.getGridVal(grid.getIndex(col - 1, row)) - currentHeight);
			heightDiff += abs(grid.getGridVal(grid.getIndex(col + 1, row)) - currentHeight);

			heightDiff += abs(grid.getGridVal(grid.getIndex(col - 1, row + 1)) - currentHeight);
			heightDiff += abs(grid.getGridVal(grid.getIndex(col, row + 1)) - currentHeight);
			heightDiff += abs(grid.getGridVal(grid.getIndex(col + 1, row + 1)) - currentHeight);

			sumHeightDiff += pow(abs(heightDiff / (currentHeight * 8)), 2);
		}
	}

	return sumHeightDiff / sqrt((double)nSamples);
}

void
CubicBezierSurfaceUtils::DumpInventorFile(std::ostream& ovstr, geometry::CubicBezierSurface &bezier, float color[3])
{
	bool convertToWorldCoordinates = true;

	if (ovstr == NULL){
		return;
	}

	Vector3D * controlPoints = bezier.GetControlPoints();

	// Write the surface out to high precision
	ovstr.precision(18);

    ovstr << "#Inventor V2.1 ascii" << std::endl;
	//  Write the header
	ovstr << "Separator {" << endl

	<< "\tComplexity {" << endl
	<< "\t\tvalue 1.0" << endl
	<< "\t}" << endl; // Complexity

	ovstr << "\tSeparator {" << endl;
	if (color != NULL){
		ovstr << "\t\tBaseColor{" << endl
		<< "\t\t\trgb " << color[0] << " " << color[1] << " " << color[2] << endl
		<< "\t\t}" << endl;
	}

	Vector3D scale (1.0, 1.0, 1.0);
	if (convertToWorldCoordinates){
		scale = geometry::Vector3D( bezier.GetPSpace().GetAxisLength( 0 ), bezier.GetPSpace().GetAxisLength( 1 ), bezier.GetPSpace().GetAxisLength( 2 ) );
	}

	ovstr << "\t\tCoordinate3 { point[" << endl;
	// patch coordinates
	for (int j = 0; j < 4; j++ ){
		for (int i = 0; i < 4; i++){
			Vector3D point = controlPoints[j * 4 + i];

			Vector3D coord = point;
			if (convertToWorldCoordinates ){
				coord = bezier.GetPSpace().getXYZfromIJK(geometry::Vector3D(point[0]/scale[0], point[1]/scale[1], point[2]/scale[2]));
			}

			ovstr << "\t\t\t" << coord[0] << " " << coord[1] << " " << coord[2] << "," << endl;
		}
	}
	// patch footer
	ovstr << "\t\t] }" << endl // Coordinate3
		<< "\t\tNurbsSurface {" << endl
		<< "\t\t\tnumUControlPoints 4" << endl
		<< "\t\t\tnumVControlPoints 4" << endl
		<< "\t\t\tuKnotVector [ 0 , 0 , 0 , 0 , 1 , 1 , 1 , 1 ]" << endl
		<< "\t\t\tvKnotVector [ 0 , 0 , 0 , 0 , 1 , 1 , 1 , 1 ]" << endl
		<< "\t\t}" << endl // NurbsSurface
		<< "\t}" << endl; // Separator 

	ovstr << "}" << endl; // Separator

	ovstr.flush();
}

void
CubicBezierSurfaceUtils::DumpPoints(std::ostream& ovstr, CubicBezierSurface & bezierCurve, double uRange[2], double vRange[2], double uInterval, double vInterval)
{
	if (ovstr == NULL){
		return;
	}

	if (uInterval <= 0.0 || uRange[0] < 0.0 || uRange[1] > 1.0 || uRange[0] > uRange[1]){
		return;
	}
	if (vInterval <= 0.0 || vRange[0] < 0.0 || vRange[1] > 1.0 || vRange[0] > vRange[1]){
		return;
	}

	ovstr.precision(18);

	for (double v = vRange[0]; v < vRange[1]; v += vInterval){
		for (double u = uRange[0]; u < uRange[1]; u += uInterval){
			Vector3D point = bezierCurve.GetPoint(u, v);
			ovstr << point[0] << " " << point[1] << " " << point[2] << endl;
		}
		Vector3D point = bezierCurve.GetPoint(uRange[1], v);
		ovstr << point[0] << " " << point[1] << " " << point[2] << endl;
	}
	for (double u = uRange[0]; u < uRange[1]; u += uInterval){
		Vector3D point = bezierCurve.GetPoint(u, vRange[1]);
		ovstr << point[0] << " " << point[1] << " " << point[2] << endl;
	}
	Vector3D point = bezierCurve.GetPoint(uRange[1], vRange[1]);
	ovstr << point[0] << " " << point[1] << " " << point[2] << endl;

	ovstr.flush();
}

void
CubicBezierSurfaceUtils::DumpInventorPoints(std::ostream& ovstr, CubicBezierSurface & bezierCurve, double uRange[2], double vRange[2], double uInterval, double vInterval, float color[3], double pointSize)
{
	bool convertToWorldCoordinates = true;

	if (ovstr == NULL){
		return;
	}

	if (uInterval <= 0.0 || uRange[0] < 0.0 || uRange[1] > 1.0 || uRange[0] > uRange[1]){
		return;
	}
	if (vInterval <= 0.0 || vRange[0] < 0.0 || vRange[1] > 1.0 || vRange[0] > vRange[1]){
		return;
	}

	ovstr.precision(18);

    ovstr << "#Inventor V2.1 ascii" << std::endl;
	//  Write the header
	ovstr << "Separator {" << endl
	<< "\tComplexity {" << endl
	<< "\t\tvalue 1.0" << endl
	<< "\t}" << endl; // Complexity

	if (color != NULL){
		ovstr << "\tBaseColor {" << endl
		<< "\t\trgb " << color[0] << " " << color[1] << " " << color[2] << " " << endl
		<< "\t}" << endl; // Complexity
	}

	ovstr << "\tSeparator {" << endl
	<< "\t\tCoordinate3 { point[" << endl;

	Vector3D scale (1.0, 1.0, 1.0);
	if (convertToWorldCoordinates){
		scale = geometry::Vector3D( bezierCurve.GetPSpace().GetAxisLength( 0 ), bezierCurve.GetPSpace().GetAxisLength( 1 ), bezierCurve.GetPSpace().GetAxisLength( 2 ) );
	}
	int count = 1;
	for (double v = vRange[0]; v < vRange[1]; v += vInterval){
		for (double u = uRange[0]; u < uRange[1]; u += uInterval){
			Vector3D point = bezierCurve.GetPoint(u, v);
			Vector3D coord = point;
			if (convertToWorldCoordinates ){
				coord = bezierCurve.GetPSpace().getXYZfromIJK(geometry::Vector3D(point[0]/scale[0], point[1]/scale[1], point[2]/scale[2]));
			}

			ovstr << "\t\t\t" << coord[0] << " " << coord[1] << " " << coord[2] << ", " << endl;

			count++;
		}
		Vector3D point = bezierCurve.GetPoint(uRange[1], v);
		Vector3D coord = point;
		if (convertToWorldCoordinates ){
			coord = bezierCurve.GetPSpace().getXYZfromIJK(geometry::Vector3D(point[0]/scale[0], point[1]/scale[1], point[2]/scale[2]));
		}

		ovstr << "\t\t\t" << coord[0] << " " << coord[1] << " " << coord[2] << ", " << endl;

		count++;
	}
	for (double u = uRange[0]; u < uRange[1]; u += uInterval){
		Vector3D point = bezierCurve.GetPoint(u, vRange[1]);
		Vector3D coord = point;
		if (convertToWorldCoordinates ){
			coord = bezierCurve.GetPSpace().getXYZfromIJK(geometry::Vector3D(point[0]/scale[0], point[1]/scale[1], point[2]/scale[2]));
		}
		ovstr << "\t\t\t" << coord[0] << " " << coord[1] << " " << coord[2] << ", " << endl;
		count++;
	}
	Vector3D point = bezierCurve.GetPoint(uRange[1], vRange[1]);
	Vector3D coord = point;
	if (convertToWorldCoordinates ){
		coord = bezierCurve.GetPSpace().getXYZfromIJK(geometry::Vector3D(point[0]/scale[0], point[1]/scale[1], point[2]/scale[2]));
	}
	ovstr << "\t\t\t" << coord[0] << " " << coord[1] << " " << coord[2] << endl
		<< "\t\t] }" << endl // Coordinate3

	<< "\t\tDrawStyle {" << endl 
	<< "\t\t\tpointSize " << pointSize << endl
	<< "\t\t}" << endl
	<< "\t\tIndexedPointSet {" << endl
	<< "\t\t\tcoordIndex [" << endl
	<< "\t\t\t\t";
	for (int i = 0; i < count; i++){
		ovstr << i << ", ";
	}
	ovstr << "-1" << endl
	<< "\t\t\t]" << endl
	<< "\t\t}" << endl // IndexedLineSet
	<< "\t}" << endl // Separator
	<< "}" << endl; // Separator
	ovstr.flush();
}


void
CubicBezierSurfaceUtils::DumpInventorControlPoints(std::ostream& ovstr, CubicBezierSurface & bezierSurface, float color[3], double pointSize)
{
	bool convertToWorldCoordinates = true;
	if (ovstr == NULL){
		return;
	}

	ovstr.precision(18);

    ovstr << "#Inventor V2.1 ascii" << std::endl;
	//  Write the header
	ovstr << "Separator {" << endl
	<< "\tComplexity {" << endl
	<< "\t\tvalue 1.0" << endl
	<< "\t}" << endl; // Complexity
	if (color != NULL){
 		ovstr << "\tBaseColor {" << endl
		<< "\t\trgb " << color[0] << " " << color[1] << " " << color[2] << " " << endl
		<< "\t}" << endl; // Complexity
	}


	ovstr << "\tSeparator {" << endl
	<< "\t\tCoordinate3 { point[" << endl;

	Vector3D * controlPoints = bezierSurface.GetControlPoints();

	Vector3D scale (1.0, 1.0, 1.0);
	if (convertToWorldCoordinates){
		scale = geometry::Vector3D( bezierSurface.GetPSpace().GetAxisLength( 0 ), bezierSurface.GetPSpace().GetAxisLength( 1 ), bezierSurface.GetPSpace().GetAxisLength( 2 ) );
	}

	for (int i = 0; i < 16; i++){
		Vector3D point = controlPoints[i];

		Vector3D coord = point;
		if (convertToWorldCoordinates ){
			coord = bezierSurface.GetPSpace().getXYZfromIJK(geometry::Vector3D(point[0]/scale[0], point[1]/scale[1], point[2]/scale[2]));
		}
		ovstr << "\t\t\t" << coord[0] << " " << coord[1] << " " << coord[2] << ", " << endl;
	}
	ovstr << "\t\t ] }" << endl;

	ovstr << "\t\tDrawStyle {" << endl 
	<< "\t\t\tpointSize " << pointSize << endl
	<< "\t\t}" << endl
	<< "\t\tIndexedPointSet {" << endl
	<< "\t\t\tcoordIndex [" << endl
	<< "\t\t\t\t";
	for (int i = 0; i < 16; i++){
		ovstr << i << ", ";
	}
	ovstr << "-1" << endl
	<< "\t\t\t]" << endl
	<< "\t\t}" << endl // IndexedLineSet
	<< "\t}" << endl // Separator
	<< "}" << endl; // Separator
	ovstr.flush();
}

void
CubicBezierSurfaceUtils::DumpInventorFile(std::ostream& ovstr, BezierQuadTree * root, float color[3], int depthLimit)
{
	bool outputControlPoints = false;

	float pointColor[3] = {1.0, 0.0, 0.0};

	if (depthLimit == -1){
		depthLimit = std::numeric_limits<int>::max();
	}

	QuadTree::QuadTreeIterator iter((QuadTree *)root);
	while ((*iter) != NULL){
		BezierQuadTree * bTree = (BezierQuadTree *)(*iter);
		if ( iter.GetDepth() > depthLimit ){
			++iter;
			continue;
		}
		if ( iter.GetDepth() != depthLimit && bTree->HaveChildren() ){
			++iter;
			continue;
		}

		CubicBezierSurface * bezierSurface = bTree->GetSurface();
			CubicBezierSurfaceUtils::DumpInventorFile(ovstr, *bezierSurface, color);

			if (outputControlPoints){
				CubicBezierSurfaceUtils::DumpInventorControlPoints(ovstr, *bezierSurface, pointColor, 4);
			}

		++iter;
	}
}

void
CubicBezierSurfaceUtils::DumpControlPointsInventorFile(std::ostream& ovstr, BezierQuadTree * root, float color[3], float pointSize, int depthLimit)
{
	if (depthLimit == -1){
		depthLimit = std::numeric_limits<int>::max();
	}

	QuadTree::QuadTreeIterator iter((QuadTree *)root);
	while ((*iter) != NULL){
		BezierQuadTree * bTree = (BezierQuadTree *)(*iter);
		if ( iter.GetDepth() > depthLimit ){
			++iter;
			continue;
		}
		if ( iter.GetDepth() != depthLimit && bTree->HaveChildren() ){
			++iter;
			continue;
		}

		CubicBezierSurface * bezierSurface = bTree->GetSurface();

		CubicBezierSurfaceUtils::DumpInventorControlPoints(ovstr, *bezierSurface, color, pointSize);

		++iter;
	}
}

void
CubicBezierSurfaceUtils::DumpPatchBoundariesInventorFile(std::ostream& ovstr, BezierQuadTree * root, float color[3], float lineThickness, int depthLimit)
{
	if (depthLimit == -1){
		depthLimit = std::numeric_limits<int>::max();
	}

	QuadTree::QuadTreeIterator iter((QuadTree *)root);
	while ((*iter) != NULL){
		BezierQuadTree * bTree = (BezierQuadTree *)(*iter);
		if ( iter.GetDepth() > depthLimit ){
			++iter;
			continue;
		}
		if ( iter.GetDepth() != depthLimit && bTree->HaveChildren() ){
			++iter;
			continue;
		}

		CubicBezierSurface * bezierSurface = bTree->GetSurface();

		Vector3D * controlPoints = bezierSurface->GetControlPoints();
		Vector3D curveControlPoints[4];
		curveControlPoints[0] = controlPoints[0];
		curveControlPoints[1] = controlPoints[1];
		curveControlPoints[2] = controlPoints[2];
		curveControlPoints[3] = controlPoints[3];
		CubicBezierCurve3D bezierCurve(curveControlPoints);
		bezierCurve.SetPSpace(bezierSurface->GetPSpace());
		CubicBezierCurve3DUtils::DumpInventorFile(ovstr, bezierCurve, color, lineThickness);

		bezierCurve.GetControlPoints()[0] = controlPoints[12 + 0];
		bezierCurve.GetControlPoints()[1] = controlPoints[12 + 1];
		bezierCurve.GetControlPoints()[2] = controlPoints[12 + 2];
		bezierCurve.GetControlPoints()[3] = controlPoints[12 + 3];
		CubicBezierCurve3DUtils::DumpInventorFile(ovstr, bezierCurve, color, lineThickness);

		bezierCurve.GetControlPoints()[0] = controlPoints[0 + 0];
		bezierCurve.GetControlPoints()[1] = controlPoints[4 + 0];
		bezierCurve.GetControlPoints()[2] = controlPoints[8 + 0];
		bezierCurve.GetControlPoints()[3] = controlPoints[12 + 0];
		CubicBezierCurve3DUtils::DumpInventorFile(ovstr, bezierCurve, color, lineThickness);

		bezierCurve.GetControlPoints()[0] = controlPoints[0 + 3];
		bezierCurve.GetControlPoints()[1] = controlPoints[4 + 3];
		bezierCurve.GetControlPoints()[2] = controlPoints[8 + 3];
		bezierCurve.GetControlPoints()[3] = controlPoints[12 + 3];
		CubicBezierCurve3DUtils::DumpInventorFile(ovstr, bezierCurve, color, lineThickness);

		++iter;
	}
}

}