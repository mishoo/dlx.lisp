;;;; package.lisp

(defpackage #:dlx
  (:use #:cl)
  (:export #:make-dlx-matrix
           #:print-dlx-solution
           #:search-dlx
           #:dlx-cover-col
           #:dlx-uncover-col
           #:dlx-node
           #:make-dlx-node
           #:dlx-node-size
           #:dlx-node-data
           #:dlx-node-up
           #:dlx-node-down
           #:dlx-node-left
           #:dlx-node-right
           #:dlx-node-col))
