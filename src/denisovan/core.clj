(ns denisovan.core
  "Namespace for the core.matrix implementation of Neanderthal.

   Should not be used directly by user code. Users should instead create and manipulate Neanderthal
   vectors and matrices using the core.matrix API"
  (:require [uncomplicate.neanderthal.core :as core] 
            [clojure.core.matrix :refer :all]
            [clojure.core.matrix.protocols :as mp]
            [uncomplicate.neanderthal.internal.host.mkl :as mkl]
            [clojure.core.matrix.implementations :as imp]
            [uncomplicate.fluokitten.core :as fluo]
            [uncomplicate.neanderthal.linalg :as lin]
            [uncomplicate.neanderthal.vect-math :as vm]
            [uncomplicate.neanderthal.math :as math]
            [uncomplicate.neanderthal.internal
             [api :as p]]
            [uncomplicate.neanderthal.vect-math :as vm])
  (:import [uncomplicate.neanderthal.internal.api Vector Matrix
            Changeable RealChangeable]))



(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(defmacro error
  "Throws an error with the provided message(s). This is a macro in order to try and ensure the 
   stack trace reports the error at the correct source line number."
  ([& messages]
    `(throw (Exception. (str ~@messages)))))

(defmacro neanderthal?
  "Returns true if v is a Neanderthal class (i.e. an instance of Matrix or Vector)"
  ([a]
    `(or (instance? Vector ~a) (instance? Matrix ~a))))

(defmacro double-coerce 
  "Coerces a 0-dimensional object to a double value"
  ([x]
  `(let [x# ~x]
     (double (if (number? x#) x# (mp/get-0d x#))))))

(defmacro tag-symbol [tag form]
  (let [tagged-sym (vary-meta (gensym "res") assoc :tag tag)]
    `(let [~tagged-sym ~form] ~tagged-sym)))

(defn vector-coerce* 
  "Coerces any numerical array to an Vector instance.
   May broadcast to the shape of an optional target if necessary.
   Does *not* guarantee a new copy - may return same data." 
  (^Vector [^Vector target m]
	  (cond (instance? Vector m) m
      (== 0 (long (mp/dimensionality m))) 
      (let [res (p/create-vector mkl/mkl-double (ecount target) false)]
        (core/transfer! (repeat (ecount target) (double-coerce m)) res))
      :else (do 
              (when (not= (ecount target) (ecount m))
                (error "Incompatible shapes coercing to vector of target length: " (ecount target)))
              (vector-coerce* m))))
  (^Vector [m]
    (cond
	    (instance? Vector m) m
      (== (dimensionality m) 1)
      (let [res (p/create-vector mkl/mkl-double (ecount m) false)]
        (core/transfer! (to-nested-vectors m) res))
      :else (error "Can't coerce to Vector with shape: " (mp/get-shape m))))) 

(defmacro vector-coerce 
  "Coerces an argument x to an Vector instance, of the same size as m"
  ([m x]
    `(tag-symbol uncomplicate.neanderthal.protocols.Vector
                 (let [x# ~x] 
                   (if (instance? Vector x#) x# (vector-coerce* ~m x#)))))
  ([x]
    `(tag-symbol uncomplicate.neanderthal.protocols.Vector
                 (let [x# ~x] 
                   (if (instance? Vector x#) x# (vector-coerce* x#))))))

(eval
  `(extend-protocol mp/PImplementation
     ~@(mapcat 
         (fn [sym]
           (cons sym
             '((implementation-key [m] :neanderthal)
               (supports-dimensionality? [m dims] (<= 1 (long dims) 2))
               (new-vector [m length] (p/create-vector mkl/mkl-double (long length) true))
               (new-matrix [m rows columns] (p/create-ge mkl/mkl-double
                                                         (long rows)
                                                         (long columns)
                                                         true
                                                         true))
                (new-matrix-nd [m shape] 
                               (case (count shape)
                                 0 0.0
                                 1 (p/create-vector mkl/mkl-double (long (first shape)) 1.0)
                                 2 (p/create-ge mkl/mkl-double
                                                (long (first shape))
                                                (long (second shape))
                                                true
                                                true)
                                 (let [d0 (first shape)
                                       moredims (next shape)]
                                   (mapv 
                                     (fn [ds] (mp/new-matrix-nd m ds))
                                     (range d0)))))
                (construct-matrix [m data]
                                  (let [dims (long (mp/dimensionality data))]
                                    (cond 
                                      (neanderthal? data) (core/copy data)
                                      (mp/is-scalar? data) (double-coerce data)
                                      (== 1 dims) (let [res (p/create-vector mkl/mkl-double (ecount data) false)]
                                                    (core/transfer! (to-nested-vectors data) res))
                                      (== 2 dims) (let [shp (shape data)
                                                        rows (long (first shp))
                                                        cols (long (second shp))
                                                        res (p/create-ge mkl/mkl-double
                                                                         rows cols
                                                                         ;; note we have to transpose since Neanderthal expects
                                                                         ;; elements in column-major order
                                                                         true
                                                                         false)]
                                                    (core/transfer! (mp/element-seq (transpose data)) res)) 
                                    :default
                                      (let [vm (mp/construct-matrix [] data)] 
                                        ;; (println m vm (shape vm))
                                         (assign! (mp/new-matrix-nd m (shape vm)) vm))))))))
         '[Vector Matrix]) ))




(extend-protocol mp/PDimensionInfo
  Vector
    (dimensionality [m]
      1)
    (is-vector? [m]
      true)
    (is-scalar? [m]
      false)
    (get-shape [m]
      [(long (.dim m))])
    (dimension-count [m x]
      (if (== 0 (long x))
        (long (.dim m))
        (error "Vector does not have dimension: " x)))
  Matrix
    (dimensionality [m]
      2)
    (is-vector? [m]
      false)
    (is-scalar? [m]
      false)
    (get-shape [m]
      [(long (.mrows m)) (long (.ncols m))])
    (dimension-count [m x]
      (let [x (long x)]
        (cond 
          (== x 0) (.mrows m)
          (== x 1) (.ncols m)
          :else (error "Matrix does not have dimension: " x)))))

(extend-protocol mp/PObjectArrayOutput
  Vector
	  (to-object-array [m]
	    (let [ec (.dim m)
	          ^objects obs (object-array ec)]
	      (dotimes [i ec] (aset obs i (.boxedEntry m i))) 
	      obs))
	  (as-object-array [m]
	    nil)
   Matrix
	  (to-object-array [m]
	    (let [rows (.mrows m)
	          cols (.ncols m)
            ^objects obs (object-array (* rows cols))]
	      (dotimes [i rows] 
          (dotimes [j cols] (aset obs (+ j (* cols i)) (.boxedEntry m i j)))) 
	      obs))
	  (as-object-array [m]
	    nil))

(extend-protocol mp/PIndexedAccess
  Vector
    (get-1d [m i]
      (.boxedEntry m (long i)))
    (get-2d [m i j]
      (error "Can't access 2-dimensional index of a vector"))
    (get-nd [m indexes]
      (when-not (== 1 (count indexes)) (error "Invalid index for Vector: " indexes))
      (.boxedEntry m (long (first indexes))))
  Matrix
    (get-1d [m i]
      (error "Can't access 1-dimensional index of a matrix"))
    (get-2d [m i j]
      (.boxedEntry m (long i) (long j)))
    (get-nd [m indexes]
      (when-not (== 2 (count indexes)) (error "Invalid index for Matrix: " indexes))
      (.boxedEntry m (long (first indexes)) (long (second indexes)))))

(extend-protocol mp/PIndexedSettingMutable
  Vector
    (set-1d! [m i v] (.set ^RealChangeable m (long i) (double-coerce v)))
    (set-2d! [m i j v] (error "Can't do 2-dimensional set on a 1D vector!"))
    (set-nd! [m indexes v]
      (if (== 1 (count indexes))
        (.set ^RealChangeable m (long (first indexes)) (double-coerce v))
        (error "Can't do " (count indexes) "-dimensional set on a 1D vector!"))) 
  Matrix
    (set-1d! [m i v] (error "Can't do 1-dimensional set on a 2D matrix!"))
    (set-2d! [m i j v] (.set ^RealChangeable m (long i) (long j) (double-coerce v)))
    (set-nd! [m indexes v]
      (if (== 2 (count indexes))
        (.set ^RealChangeable m (long (first indexes)) (long (second indexes)) (double-coerce v))
        (error "Can't do " (count indexes) "-dimensional set on a 2D matrix!"))))

(extend-protocol mp/PSubVector
  Vector
    (subvector [m start length]
      (let [k (long start)
            l (long length)]
        (.subvector m k l)))) 

(extend-protocol mp/PSliceView
  Matrix
    (get-major-slice-view [m i] 
      (.row m (long i))))

(extend-protocol mp/PSliceView2
  Matrix
    (get-slice-view [m dim i]
      (case (long dim) 
        0 (.row m (long i))
        1 (.col m (long i))
        (error "Can't slice on dimension " dim " of a Matrix"))))

(extend-protocol mp/PSliceSeq
  Matrix  
    (get-major-slice-seq [m] 
      (mapv (fn [i] (.row m (long i)))
            (range (.mrows m)))))

(extend-protocol mp/PMatrixRows
  Matrix
    (get-rows [m]
      (mapv (fn [i] (.row m (long i)))
            (range (.mrows m)))))

(extend-protocol mp/PMatrixColumns
  Matrix
    (get-columns [m]
      (mapv (fn [i] (.col m (long i)))
            (range (.ncols m)))))

(extend-protocol mp/PTranspose
  Vector (transpose [m] m)
  Matrix (transpose [m] (.transpose m))) 

(extend-protocol mp/PMatrixCloning
  Vector (clone [m] (core/copy m))
  Matrix	(clone [m] (core/copy m)))

(extend-protocol mp/PElementCount
  Matrix
    (element-count [m]
      (.dim m))
  Vector
    (element-count [m]
      (.dim m)))

(extend-protocol mp/PSummable
  Vector
    (element-sum [m]
      (core/sum m))
  Matrix
    (element-sum [m]
      (reduce + (map core/sum (slices m)))))

(extend-protocol mp/PVectorOps
  Vector
    (vector-dot [a b]
      (core/dot a (vector-coerce a b)))
    (length [a]
      (core/nrm2 a))
    (length-squared [a]
      (let [l ^double (core/nrm2 a)]
        (* l l)))
    (normalise [a]
      (let [l (core/nrm2 a)]
        (div a l))))

;; fixes test
(extend-protocol mp/PSquare
  Vector
  (square [m] (vm/sqr m))
  Matrix
  (square [m] (vm/sqr m)))

;; fixes test

(eval
 `(extend-protocol mp/PMatrixDivide
    ~@(mapcat
        (fn [sym]
          (cons
           sym
           '((element-divide [m] (vm/inv m))
             (element-divide [m a]
                             (if (number? a)
                               (let [a (double a)]
                                 (core/scal (/ 1.0 a) m))
                               (let [[m a] (mp/broadcast-compatible m a)]
                                 (vm/div m a)))))))
        '[Vector Matrix])))

;; optional, for performance

;; make it fast
(extend-protocol mp/PMatrixMultiply
  Matrix
  (matrix-multiply [m a]
    (cond (matrix? a)
          (core/mm m (matrix a))
          (vec? a)
          (core/mv m (matrix a))))
  (element-multiply [m a]
    (if (number? a)
      (let [a (double a)]
        (core/scal a m))
      ;; TODO broadcasting perf. critical
      (let [[m a] (mp/broadcast-compatible m a)]
        (vm/mul m a)))))


(extend-protocol mp/PMatrixProducts
  Vector
  (inner-product [m a]
    (if (vec? a)
      (core/dot m a)
      (core/mv a m)))
  (outer-product [m a]
    (when (vec? a)
      (core/rk m a))
    ;; else fall-back
    ))

(eval
 `(extend-protocol mp/PAddScaled
    ~@(mapcat
        (fn [sym]
          (cons
           sym
           '((add-scaled [m a factor]
               (core/axpy factor a m)))))
        '[Vector Matrix])))

(eval
 `(extend-protocol mp/PAddScaledMutable
    ~@(mapcat
        (fn [sym]
          (cons
           sym
           '((add-scaled! [m a factor]
                         (core/axpy! factor a m)))))
        '[Vector Matrix])))


(extend-protocol mp/PAddScaledMutable
  Vector
  (add-scaled [m a factor]
    (core/axpy factor a m)))


(defn svd-inv [m]
  (let [{:keys [u vt sigma] :as svd} (lin/svd m true true)
        diag-inverse (core/gd mkl/mkl-double (core/mrows sigma)
                              ;; TODO seq somehow  needed?
                              (seq (fluo/fmap (fn ^double [^double x] (/ 1 x)) sigma)))]
    (core/mm (core/trans vt) diag-inverse (core/trans u))))


(extend-protocol mp/PMatrixOps
  Matrix
  (trace [m] (core/sum (core/dia m)))
  (det [m] (lin/det m))
  (inverse [m] (svd-inv m)))


(eval
 `(extend-protocol mp/PExponent
    ~@(mapcat
        (fn [sym]
          (cons
           sym
           '((element-pow [m exponent]
                          (vm/pow m exponent)))))
        '[Vector Matrix])))

(eval
 `(extend-protocol mp/PLogistic
    ~@(mapcat
       (fn [sym]
         (cons
          sym
          '((logistic [m]
                      (let [e (vm/exp (core/scal -1.0 m))]
                       (vm/inv (vm/linear-frac e 1.0)))))))
       '[Vector Matrix])))

(eval
 `(extend-protocol mp/PMathsFunctions
    ~@(mapcat
        (fn [sym]
          (cons
           sym
           '((abs [m] (vm/abs m))
             (acos [m] (vm/acos m))
             (asin [m] (vm/asin m))
             (atan [m] (vm/atan m))
             (cbrt [m] (vm/cbrt m))
             (ceil [m] (vm/ceil m))
             (cos [m] (vm/cos m))
             (cosh [m] (vm/cosh m))
             (exp [m] (vm/exp m))
             (floor [m] (vm/floor m))
             (log [m] (vm/log m))
             (log10 [m] (vm/log10 m))
             (round [m] (vm/round m))
             #_(signum [m] (vm/si))
             (sin [m] (vm/sin m))
             (sinh [m] (vm/sinh m))
             (tan [m] (vm/tan m))
             (tanh [m] (vm/tanh m))
             #_(to-degrees [m] (vm/to))
             #_(to-radians [m]))))
        '[Vector Matrix])))

(eval
 `(extend-protocol mp/PMathsFunctionsMutable
    ~@(mapcat
        (fn [sym]
          (cons
           sym
           '((abs! [m] (vm/abs! m))
             (acos! [m] (vm/acos! m))
             (asin! [m] (vm/asin! m))
             (atan! [m] (vm/atan! m))
             (cbrt! [m] (vm/cbrt! m))
             (ceil! [m] (vm/ceil! m))
             (cos! [m] (vm/cos! m))
             (cosh! [m] (vm/cosh! m))
             (exp! [m] (vm/exp! m))
             (floor! [m] (vm/floor! m))
             (log! [m] (vm/log! m))
             (log10! [m] (vm/log10! m))
             (round! [m] (vm/round! m))
             #_(signum [m] (vm/si))
             (sin! [m] (vm/sin! m))
             (sinh! [m] (vm/sinh! m))
             (tan! [m] (vm/tan! m))
             (tanh! [m] (vm/tanh! m))
             #_(to-degrees [m] (vm/to))
             #_(to-radians [m]))))
        '[Vector Matrix])))


;; Register the Neanderthal implementation using MKL
(imp/register-implementation (p/create-vector mkl/mkl-double 3 true))


(comment
  (set-current-implementation :neanderthal)
  (set-current-implementation :persistent-vector)

  (add-scaled (matrix [[1 2] [3 4]]) (matrix [[1 2] [3 4]]) 0.1)

  (add-scaled (matrix [[1 2] [3 4]]) (matrix [[1 2] [3 4]]) 0.1)

  (outer-product (matrix [[1 2] [3 4]]) (matrix [[1 2] [3 4]]))


  (core/ge (core/vctr [1 2 3]))

  (vector-coerce [1 2 3])

  (mmul (identity-matrix 3) [1 2 3])

  (require '[uncomplicate.neanderthal.linalg :as lin])

  (broadcast (matrix [1 2 3]) [10 3])

  (core/copy! (matrix [1 2 3]) (matrix [0 0 0 0 0 0]) 1 2 2)

  (core/vctr mkl/mkl-double (matrix [[1 2 3]]))

  (core/copy! (matrix [[1 2 3]])
              (matrix [[0 0 0 0]]) 0 3 0)



  (let [test-mtx (matrix [1 2 3]) 
        I (identity-matrix 3)]
    (equals test-mtx (mmul I test-mtx))
    )

  

  (let [m (matrix [[1 2 3]
                   [2 3 4]])]
    (div  m  [1 2 3])
    #_(equals m  1.0E-4))

  (require '[clojure.core.matrix.operators :as ops])

  (let [m (matrix [[1 2 3]])]
    (equals (square m) (ops/** m 2) 1.0E-4))

  (div (matrix [[1 0] [0 1]]) 3)


  (let [A (identity-matrix 1000)
        B (identity-matrix 1000)]
    (time (do (mmul A B)
              #_(core/mm A B)
              :done)))

  (let [a (matrix (range 100000))
        b (matrix (reverse (range 100000)))]
    (time
     (dot a b)
     #_(core/dot a b)))

  (mmul (matrix [[1 0] [0 1]]) (matrix [1 0]) )

  

  (exp (matrix [1 2 3])))



