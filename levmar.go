package problems

/*
#cgo LDFLAGS: -L/usr/lib  -L/home/tony/src/levmar-2.6 -llevmar -llapack -lblas -lf2c  -lm

void func(double *p, double *x, int m, int n, void *data);
void jacfunc(double *p, double *x, int m, int n, void *data);
void levmar(double* ygiven, double* p, const int n, const int m, void* e );

*/
import "C"

import (
	expr "damd/go-symexpr"
	"math"
	"reflect"
	"unsafe"
)

type callback_data struct {
	Train []*PointSet
	Test  []*PointSet
	E     expr.Expr
	J     []expr.Expr
	Coeff []float64
	Task  ExprProblemType
}

func RegressExpr(E expr.Expr, P *ExprProblem) (R *ExprReport) {

	c := make([]float64, 0)
	c, eqn := E.ConvertToConstants(c)

	var coeff []float64
	if len(c) > 0 {
		coeff = LevmarExpr(eqn, P.SearchVar, P.SearchType, c, P.Train, P.Test)
	}

	R = new(ExprReport)
	R.expr = eqn /*.ConvertToConstantFs(coeff)*/
	R.coeff = coeff
	_, s2, serr := scoreExpr(E, P, coeff)
	R.testScore = s2
	R.testError = serr
	R.expr.CalcExprStats(0)

	return R
}

func scoreExpr(e expr.Expr, P *ExprProblem, coeff []float64) (int, int, float64) {
	score := 0
	score2 := 0
	error := 0.0

	for _, PS := range P.Test {
		for _, p := range PS.Points() {
			y := p.Depnd(P.SearchVar)
			var out float64
			if P.SearchType == ExprBenchmark {
				out = e.Eval(0, p.Indeps(), coeff, PS.SysVals())
			} else if P.SearchType == ExprDiffeq {
				out = e.Eval(p.Indep(0), p.Indeps()[1:], coeff, PS.SysVals())
			}

			diff := math.Abs(out - y)
			if math.IsNaN(diff) {
				continue
			}
			if diff < P.HitRatio {
				score++
			}
			err := math.Abs(diff / y)
			if math.IsNaN(err) || math.IsInf(err, 0) {
				err = diff
			}
			error += err
			if err < P.HitRatio {
				score2++
			}
		}
	}

	eAve := error / (float64(len(P.Test)) * float64(P.Test[0].NumPoints()))
	// eAve := error / float64(P.Test.NumPoints())

	return score, score2, eAve
}

func LevmarExpr(e expr.Expr, searchDim int, task ExprProblemType, guess []float64, train, test []*PointSet) []float64 {

	ps := train[0].NumPoints()
	PS := len(train) * ps

	c := make([]float64, len(guess))
	var cd callback_data
	cd.Train = train
	cd.Test = test
	cd.E = e
	cd.Coeff = c
	cd.Task = task
	cd.J = make([]expr.Expr, len(guess))
	for i, _ := range cd.J {
		cd.J[i] = e.DerivConst(i) /*.Simplify(SimpRules{true,true})*/
	}

	// c/levmar inputs
	coeff := make([]C.double, len(guess))
	for i, g := range guess {
		coeff[i] = C.double(g)
	}

	y := make([]C.double, PS)
	for i1, T := range train {
		for i2, p := range T.Points() {
			i := i1*ps + i2
			y[i] = C.double(p.Depnd(searchDim))
		}
	}
	ya := (*C.double)(unsafe.Pointer(&y[0]))
	ca := (*C.double)(unsafe.Pointer(&coeff[0]))
	ni := C.int(PS)
	mi := C.int(len(c))

	C.levmar(ya, ca, ni, mi, unsafe.Pointer(&cd))

	for i, _ := range coeff {
		c[i] = float64(coeff[i])
	}
	return c
}

//export Callback_func
func Callback_func(p, x *C.double, e unsafe.Pointer) {

	cd := *(*callback_data)(e)
	eqn := cd.E
	coeff := cd.Coeff

	M1 := len(cd.Train)
	M2 := cd.Train[0].NumPoints()
	M3 := len(coeff)
	M23 := M2 * M3
	MA := M1 * M23

	var p_go []C.double
	p_head := (*reflect.SliceHeader)((unsafe.Pointer(&p_go)))
	p_head.Cap = M3
	p_head.Len = M3
	p_head.Data = uintptr(unsafe.Pointer(p))
	for i, _ := range p_go {
		coeff[i] = float64(p_go[i])
	}

	var x_go []C.double
	x_head := (*reflect.SliceHeader)((unsafe.Pointer(&x_go)))
	x_head.Cap = MA
	x_head.Len = MA
	x_head.Data = uintptr(unsafe.Pointer(x))

	N := len(cd.Train)
	var out float64
	for i1, PS := range cd.Train {
		for i2, pnt := range PS.Points() {
			i := i1*N + i2
			if cd.Task == ExprBenchmark {
				out = eqn.Eval(0, pnt.Indeps(), coeff, PS.SysVals())
			} else if cd.Task == ExprDiffeq {
				out = eqn.Eval(pnt.Indep(0), pnt.Indeps()[1:], coeff, PS.SysVals())
			}

			x_go[i] = C.double(out)

		}
	}
}

//export Callback_jacfunc
func Callback_jacfunc(p, x *C.double, e unsafe.Pointer) {

	cd := *(*callback_data)(e)
	coeff := cd.Coeff

	M1 := len(cd.Train)
	M2 := cd.Train[0].NumPoints()
	M3 := len(coeff)
	M23 := M2 * M3
	MA := M1 * M23

	var p_go []C.double
	p_head := (*reflect.SliceHeader)((unsafe.Pointer(&p_go)))
	p_head.Cap = M3
	p_head.Len = M3
	p_head.Data = uintptr(unsafe.Pointer(p))
	for i, _ := range p_go {
		coeff[i] = float64(p_go[i])
	}

	var x_go []C.double
	x_head := (*reflect.SliceHeader)((unsafe.Pointer(&x_go)))
	x_head.Cap = MA
	x_head.Len = MA
	x_head.Data = uintptr(unsafe.Pointer(x))

	var out float64
	for i1, PS := range cd.Train {
		for i2, pnt := range PS.Points() {
			i := i1*M23 + i2*M3
			for ji, eqn := range cd.J {
				if cd.Task == ExprBenchmark {
					out = eqn.Eval(0, pnt.Indeps(), coeff, PS.SysVals())
				} else if cd.Task == ExprDiffeq {
					out = eqn.Eval(pnt.Indep(0), pnt.Indeps()[1:], coeff, PS.SysVals())
				}

				x_go[i+ji] = C.double(out)
			}
		}
	}

}
