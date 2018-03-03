(ns denisovan.core-test
  (:use [clojure test]
        [clojure.core.matrix])
  (:require [clojure.core.matrix.protocols :refer [norm] :as p]
            [denisovan.core :as d]
            [clojure.core.matrix.compliance-tester]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(set-current-implementation :neanderthal)

(deftest instance-tests
  (clojure.core.matrix.compliance-tester/instance-test (array :neanderthal [[1 2 3] [4 5 6]]))
  (clojure.core.matrix.compliance-tester/instance-test (array :neanderthal [1 2 3 4 5]))) 

(deftest compliance-test
  (clojure.core.matrix.compliance-tester/compliance-test (array :neanderthal [1 2]))) 

(def eps 1e-8)

(deftest additional-operations
  (testing "Operations tested additional to compliance tests."
    (is (<= (norm (sub (sub 1.0 (matrix [1 2 3]))
                      (matrix [0 -1 -2]))
                 2)
           eps))
    (is (<= (norm (sub (scale-add! (matrix [1 2]) 3 (matrix [11 17]) 5 0.0)
                      (matrix [58 91]))
                  2)
            eps))

    (let [m (matrix [[1 2] [3 4]])]
      (is (<= (norm (sub (mmul (d/svd-inv m) m)
                       (identity-matrix 2))
                   2)
             eps)))))
