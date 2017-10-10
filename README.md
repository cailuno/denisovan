# denisovan

This library provides a [core.matrix](https://github.com/mikera/core.matrix)
implementation for [neanderthal](http://neanderthal.uncomplicate.org/). The main
focus of this library is to map neanderthal's high performance BLAS routines to
core.matrix protocols as closely as possible, while being compliant with the
rest of core.matrix. 

Please take a look at neanderthal if you want to improve your performance
further, as it is exposing the high performance low-level primitives directly in
Clojure and provides many knobs to implement fast numerical algorithms from the
literature. If you just need standard matrix multiplications and operations, as
e.g. in deep-learning or general optimization algorithms you should be able to
use denisovan just fine without losing portability to other core.matrix
backends. If in doubt use criterium and a direct implementation with neanderthal
to check whether your operations do something inefficiently behind the
core.matrix API. Please join
our [gitter](https://gitter.im/metasoarous/clojure-datascience) chat on
questions and feel free to open issues!

## Usage

You just need to load and activate the implementation as usual.

~~~clojure
(require '[denisovan.core])

(clojure.core.matrix/set-current-implementation :neanderthal)
~~~

## License

Copyright © 2016 Mike Anderson, 2017 Christian Weilbach

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.
