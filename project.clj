(defproject denisovan "0.1.0-SNAPSHOT"
  :description "A core.matrix backend for neanderthal."
  :url "http://github.com/whilo/denisovan"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.9.0"]
                 [uncomplicate/neanderthal "0.18.0"]
                 [net.mikera/core.matrix "0.62.0"]
                 [uncomplicate/fluokitten "0.6.1"]
                 #_[plotly-clj "0.1.1"]]

  :profiles {:dev {:global-vars {*warn-on-reflection* true
                                 *assert* false
                                 *unchecked-math* :warn-on-boxed
                                 *print-length* 128}
                   :dependencies [[net.mikera/core.matrix "0.62.0" :classifier "tests"]
                                  [net.mikera/vectorz-clj "0.47.0" :exclusions [net.mikera/core.matrix]]
                                  [criterium/criterium "0.4.4"]]}}
  )
