package problems

/*
#cgo LDFLAGS: -L/usr/lib  -Llevmar-2.6 -llevmar -llapack -lblas -lf2c  -lm

#include "levmar_h.h"
#include "stack.h"

*/
import "C"

import (
	// "fmt"
	expr "github.com/verdverm/go-symexpr"
	. "pge1/problems"
	"reflect"
	"unsafe"
)

type C_double C.double

func MakeCDouble(f float64) C_double {
	return C_double(C.double(f))
}

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

	// C.levmar_dif(ya, ca, mi, ni, unsafe.Pointer(&cd))
	C.levmar_der(ya, ca, mi, ni, unsafe.Pointer(&cd))

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

func StackLevmarExpr(e expr.Expr, x_dims int, coeff []float64, c_ygiven, c_input []C_double) []float64 {

	// fmt.Printf("Entering Stack Levmar: \n")

	c_coeff := make([]C.double, len(coeff))
	for i, c := range coeff {
		c_coeff[i] = C.double(c)
	}
	// c_ygiven := make([]C.double, len(ygiven))
	// for i, c := range ygiven {
	// 	c_ygiven[i] = C.double(c)
	// }
	// c_input := make([]C.double, len(input))
	// for i, c := range input {
	// 	c_input[i] = C.double(c)
	// }

	cp := (*C.double)(unsafe.Pointer(&c_coeff[0]))
	mi := C.int(len(coeff))
	yp := (*C.double)(unsafe.Pointer(&c_ygiven[0]))
	ni := C.int(len(c_ygiven))

	var sd *C.StackData
	sd = new(C.StackData)
	// fmt.Printf("x_len: %d   x_dim: %d\n", len(input), x_dims)
	sd.x_len = C.int(len(c_input))
	sd.x_dim = C.int(x_dims)
	sd.x_data = (*C.double)(unsafe.Pointer(&c_input[0]))

	serial := make([]int, 0)
	serial = e.StackSerial(serial)
	// fmt.Printf("StackSerial: %v\n", serial)
	c_serial := make([]C.int, len(serial))
	for i, I := range serial {
		c_serial[i] = C.int(I)
	}
	sd.expr.serial = (*C.int)(unsafe.Pointer(&c_serial[0]))
	sd.expr.s_len = C.int(len(serial))

	// fmt.Printf("GOT HERE: %v\n", serial)

	derivs := make([]C.StackExpr, len(coeff))
	for i, _ := range coeff {
		deriv := e.DerivConst(i)
		dserial := make([]int, 0)
		dserial = deriv.StackSerial(dserial)
		d_serial := make([]C.int, len(dserial))
		for i, I := range dserial {
			d_serial[i] = C.int(I)
		}
		derivs[i].serial = (*C.int)(unsafe.Pointer(&d_serial[0]))
		derivs[i].s_len = C.int(len(d_serial))

	}
	sd.derivs = (*C.StackExpr)(unsafe.Pointer(&derivs[0]))
	sd.d_len = C.int(mi)

	// fmt.Printf("Going into stack_levmar_dif\n")
	C.stack_levmar_dif(yp, cp, mi, ni, unsafe.Pointer(sd))
	// C.stack_levmar_der(yp, cp, mi, ni, unsafe.Pointer(sd))
	// fmt.Printf("Returned from stack_levmar_dif\n")

	c := make([]float64, len(c_coeff))
	for i, _ := range c_coeff {
		c[i] = float64(c_coeff[i])
	}
	// fmt.Printf("C0: %f\n", c[0])
	return c
}

// ./pge1 -pcfg=prob/bench/Koza_1.cfg -peel=3 -iter=100 -init=method1 -grow=method1
