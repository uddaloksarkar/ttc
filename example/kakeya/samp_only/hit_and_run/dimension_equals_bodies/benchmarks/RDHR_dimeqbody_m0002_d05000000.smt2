(set-logic QF_LRA)
(declare-fun x0 () Real)
(declare-fun x1 () Real)
(assert (or
  (and
    (<= (- 1) x0)
    (<= x0 1)
    (<= (- 100) x1)
    (<= x1 100)
  )
  (and
    (<= (- 100) x0)
    (<= x0 100)
    (<= (- 1) x1)
    (<= x1 1)
  )
))
(check-sat)
