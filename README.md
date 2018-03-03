# denisovan

This library provides a [core.matrix](https://github.com/mikera/core.matrix)
implementation for [neanderthal](http://neanderthal.uncomplicate.org/). The main
focus of this library is to map neanderthal's high performance BLAS routines to
core.matrix protocols as closely as possible, while being compliant with the
rest of core.matrix. For documentation please look at the `core.matrix` and
neanderthal documentation, as this library is mostly glue code. If you encounter
any issues, including surprising performance behaviour, please open an issue.

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

Include this in your dependencies:

[![Clojars Project](https://img.shields.io/clojars/v/denisovan.svg)](https://clojars.org/denisovan)

You just need to load and activate the implementation as usual.

~~~clojure
(require '[denisovan.core])

(clojure.core.matrix/set-current-implementation :neanderthal)
~~~

## Roadmap

### 0.2.0
- support GPU and OpenCL backend
- add any additional protocols that are important for performance

## License

Copyright © 2016 Mike Anderson, 2017 Christian Weilbach

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.
