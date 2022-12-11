;;;; This package implements Knuth's Algorithm X using the Dancing
;;;; Links technique, which he suggests we call the DLX algorithm.
;;;;
;;;; https://www-cs-faculty.stanford.edu/~uno/papers/dancing-color.ps.gz
;;;;
;;;; You should read the paper (at least the first 6 pages) to get
;;;; familiar with the ideas, in order for this package to be of any
;;;; use.

(in-package #:dlx)

(defstruct dlx-node
  "A node in Knuth's fabulous data structure. For simplicity we use the
same structure for all nodes (including headers), even though,
depending on the node type, some fields might be unused."
  (size 0 :type fixnum)
  (data nil)
  (up nil :type (or null dlx-node))
  (down nil :type (or null dlx-node))
  (left nil :type (or null dlx-node))
  (right nil :type (or null dlx-node))
  (col nil :type (or null dlx-node)))

(defun append-horiz (prev new)
  "Inserts `new' node horizontally in the circular doubly-linked list,
after `prev'."
  (setf (dlx-node-right new) (dlx-node-right prev)
        (dlx-node-left (dlx-node-right prev)) new
        (dlx-node-right prev) new
        (dlx-node-left new) prev))

(defun insert-vert (anchor new)
  "Inserts `new' node vertically in the circular doubly-linked list,
before `anchor'."
  (setf (dlx-node-down new) anchor
        (dlx-node-up new) (dlx-node-up anchor)
        (dlx-node-down (dlx-node-up anchor)) new
        (dlx-node-up anchor) new))

(defun make-dlx-matrix (matrix &optional coldata)
  "Creates Knuth's toroidal data structure, given `matrix' of 1s and 0s,
and optional `coldata' (must be a list if present) for the column
headers `data' slots. `matrix' must be a vector of rows, where rows
are vectors of 1 and 0 (each row could be a bit vector). All rows
should be of the same length. Example #(#*101 #*010).

Returns `head', the matrix entry point (the root object, as Knuth
calls it, which is noted `h' in his paper)."
  (let ((ncols (length (aref matrix 0)))
        (head (make-dlx-node)))
    (loop repeat ncols
          for data = coldata then (cdr data)
          for prev = head then node
          for node = (make-dlx-node :data (car data))
          do (setf (dlx-node-right prev) node
                   (dlx-node-left node) prev
                   (dlx-node-up node) node
                   (dlx-node-down node) node)
          finally
             (setf (dlx-node-left head) node
                   (dlx-node-right node) head))
    (loop for i from 0
          for row across matrix
          do (loop with new and prev = nil
                   for col = (dlx-node-right head) then (dlx-node-right col)
                   for bit across row
                   unless (zerop bit)
                     do (setf new (make-dlx-node :col col :data i))
                        (if prev
                            (append-horiz prev new)
                            (setf (dlx-node-left new) new
                                  (dlx-node-right new) new))
                        (insert-vert col new)
                        (incf (dlx-node-size col))
                        (setf prev new)))
    head))

(defmacro find-best (data predicate exp &body init)
  (setf predicate (case predicate
                    (minimizing '<)
                    (maximizing '>)
                    (otherwise predicate)))
  (let ((best-data (gensym "best-data"))
        (best-val (gensym "best-val"))
        (current-data (gensym "current-data"))
        (current-val (gensym "current-val")))
    `(loop with ,best-data = nil and ,best-val = nil ,@init
           for ,current-data = ,data
           for ,current-val = ,exp
           when (or (not ,best-val)
                    (,predicate ,current-val ,best-val))
             do (setf ,best-val ,current-val
                      ,best-data ,current-data)
           finally (return (values ,best-data ,best-val)))))

(defun select-col (head)
  ;; wish `loop' could do this..
  (find-best col minimizing (dlx-node-size col)
    for col = (dlx-node-right head) then (dlx-node-right col)
    until (eq col head)))

(defun print-dlx-solution (sol)
  "The default solution printer. Returns the solution selected by the
algorithm as a list of indexes of the rows that fit the constraints."
  (mapcar #'dlx-node-data sol))

(defun search-dlx (head &key
                          (solcount nil)
                          (printer 'print-dlx-solution))
  "Search for solutions. `head' is the root object. Pass
`solcount' (integer) if you want to limit the number of solutions (by
default it searches for all solutions, which is impractical for many
problems). `printer' should be `nil', or a function that converts the
solution. The default printer will convert to row indexes. The print
function receives a list of nodes selected by the algorithm. Check
Knuth's paper to see how he suggests printing solution (page 5 at the
bottom).

If you pass `nil' for the `printer', this function will only count the
solutions, which is faster and should not do any consing. Otherwise it
returns a list of solutions.

Unless interrupted, the matrix will not be altered; even when the
search is limited by `solcount', links will keep dancing until they
get back to the original shape."
  (let ((solutions nil)
        (found 0))
    (labels ((run (sol)
               (when (eq (dlx-node-right head) head)
                 (when printer
                   (push (funcall printer sol) solutions))
                 (incf found)
                 (return-from run))
               (let ((c (select-col head)))
                 (dlx-cover-col c)
                 (loop for r = (dlx-node-down c) then (dlx-node-down r)
                       until (eq r c)
                       until (and solcount (= found solcount))
                       do (loop for j = (dlx-node-right r) then (dlx-node-right j)
                                until (eq j r)
                                do (dlx-cover-col (dlx-node-col j)))
                          (run (when printer (cons r sol)))
                          (loop for j = (dlx-node-left r) then (dlx-node-left j)
                                until (eq j r)
                                do (dlx-uncover-col (dlx-node-col j))))
                 (dlx-uncover-col c))))
      (run nil))
    (if printer solutions found)))

(defun dlx-cover-col (c)
  "Covers the given column."
  (setf (dlx-node-left (dlx-node-right c)) (dlx-node-left c)
        (dlx-node-right (dlx-node-left c)) (dlx-node-right c))
  (loop for i = (dlx-node-down c) then (dlx-node-down i)
        until (eq i c)
        do (loop for j = (dlx-node-right i) then (dlx-node-right j)
                 until (eq j i)
                 do (setf (dlx-node-up (dlx-node-down j)) (dlx-node-up j)
                          (dlx-node-down (dlx-node-up j)) (dlx-node-down j))
                    (decf (dlx-node-size (dlx-node-col j))))))

(defun dlx-uncover-col (c)
  "Uncovers a column that was previously covered."
  (loop for i = (dlx-node-up c) then (dlx-node-up i)
        until (eq i c)
        do (loop for j = (dlx-node-left i) then (dlx-node-left j)
                 until (eq j i)
                 do (incf (dlx-node-size (dlx-node-col j)))
                    (setf (dlx-node-up (dlx-node-down j)) j
                          (dlx-node-down (dlx-node-up j)) j)))
  (setf (dlx-node-left (dlx-node-right c)) c
        (dlx-node-right (dlx-node-left c)) c))
