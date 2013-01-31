package problems

/*
#cgo LDFLAGS: -L/usr/lib  -Llevmar-2.6 -llevmar -llapack -lblas -lf2c  -lm

#include "levmar_h.h"

*/
import "C"

import (
	expr "github.com/verdverm/go-symexpr"
	. "pge1/problems"
	"reflect"
	"unsafe"
)

type callback_data struct {
	Train []*PointSet
	Test  []*PointSet
	E     expr.Expr
	J     []expr.Expr
	S     [][]int
	Coeff []float64
	Task  ExprProblemType
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
		deriv := e.DerivConst(i)
		// simp := deriv.Simplify(expr.DefaultRules())
		cd.J[i] = deriv
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


	C.levmar(ya, ca, mi, ni, unsafe.Pointer(&cd))

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




func LevmarStack(e expr.Expr, searchDim int, task ExprProblemType, guess []float64, train, test []*PointSet) []float64 {

	ps := train[0].NumPoints()
	PS := len(train) * ps
	x_dim := train[0].NumDim()
	x_len := PS
	x_tot := PS * x_dim

	// c/levmar inputs
	coeff := make([]C.double, len(guess))
	for i, g := range guess {
		coeff[i] = C.double(g)
	}

	x := make([]C.double, x_tot)
	y := make([]C.double, PS)
	for i1, T := range train {
		for i2, p := range T.Points() {
			i := i1*ps + i2
			y[i] = C.double(p.Depnd(searchDim))
			for i3, x_p := range p.Indeps() {
				j := i1*x_len + i2*T.NumPoints() + i3
				x[j] = C.double(x_p)
			}
		}
	}
	ya := (*C.double)(unsafe.Pointer(&y[0]))
	ca := (*C.double)(unsafe.Pointer(&coeff[0]))
	ni := C.int(PS)
	mi := C.int(len(guess))

	var sd *C.StackData
	sd.x_len = C.int(x_tot)
	sd.x_dim = C.int(x_dim)
	sd.d_len = C.int(mi)
	sd.x_data = (*C.double)(unsafe.Pointer(&x[0]))


	serial := make([]int,0)
	serial = e.StackSerial(serial)
	sd.expr.s_len = C.int(len(serial))
	c_serial := make([]C.int,len(serial))
	for i, I := range serial {
		c_serial[i] = C.int(I)
	}
	sd.expr.serial = (*C.int)(unsafe.Pointer(&c_serial[0]))

	derivs := make([]C.StackExpr,len(guess))
	for i, _ := range guess {
		deriv := e.DerivConst(i)
		serial := make([]int,0)
		serial = deriv.StackSerial(serial)
		derivs[i].s_len = C.int(len(serial))
		c_serial := make([]C.int,len(serial))
		for i, I := range serial {
			c_serial[i] = C.int(I)
		}
		derivs[i].serial = (*C.int)(unsafe.Pointer(&c_serial[0]))

	}
	sd.derivs = (*C.StackExpr)(unsafe.Pointer(&derivs[0]))


	C.stack_levmar(ya, ca, mi, ni, unsafe.Pointer(&sd))


	c := make([]float64, len(guess))
	for i, _ := range coeff {
		c[i] = float64(coeff[i])
	}
	return c
}


