;;;; A Clojure implementation of the Mersenne Twister algorithm  
;;;;
;;;; Copyright (c) 2009 Steve Knight <stkni@gmail.com>
;;;;
;;;; Copyright (C) 1997, 1999 Makoto Matsumoto and Takuji Nishimura. 
;;;; matumoto@math.keio.ac.jp    
;;;; 
;;;; Permission is hereby granted, free of charge, to any person obtaining
;;;; a copy of this software and associated documentation files (the
;;;; "Software"), to deal in the Software without restriction, including
;;;; without limitation the rights to use, copy, modify, merge, publish,
;;;; distribute, sublicense, and/or sell copies of the Software, and to
;;;; permit persons to whom the Software is furnished to do so, subject to
;;;; the following conditions:
;;;; 
;;;; The above copyright notice and this permission notice shall be
;;;; included in all copies or substantial portions of the Software.
;;;; 
;;;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
;;;; EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
;;;; MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
;;;; NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
;;;; LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
;;;; OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
;;;; WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
;;;;
                                    
(ns com.hackinghat.common-library.mersenne-twister)

(def n 624)
(def m 397)
(def matrix-a 0x9908b0df)            ; constant vector a
(def upper-mask 0x80000000)           ; most significant w-r bits
(def lower-mask 0x7fffffff)           ; least significant w-r bits
(def tempering-mask-b 0x9d2c5680)
(def tempering-mask-c 0xefc60000)
(defn tempering-shift-u [y] (bit-shift-right y 11))
(defn tempering-shift-s [y] (bit-shift-left y 7))
(defn tempering-shift-t [y] (bit-shift-left y 15))
(defn tempering-shift-l [y] (bit-shift-right y 18))

(def mt  (ref []))
(def mti (ref (inc n)))

(defn reset []
  (dosync 
    (ref-set mt [])
    (ref-set mti (inc n)))
  nil)
    
(defn sgenrand [seed]
  "Seed the original array"
  (dosync
    (ref-set mt [])
    (loop [i (int 0) s seed]
      (if (< i n)
        (let [mti     (bit-and s 0xffff0000)
              s-1     (inc (* 69069 s))]
           (ref-set mt (into @mt [(bit-or mti (bit-shift-right (bit-and s-1 0xffff0000) 16))]))
           (recur (inc i) (inc (* 69069 s-1))))))
    (ref-set mti n)))

(defn -genrand-inner [v1 v2 v3]
  "Factored from the reference implementation"
  (let [y (bit-or (bit-and v1 upper-mask)
                  (bit-and v2 lower-mask))]
    (bit-xor (bit-xor v3 (bit-shift-right y 1))
             (if (= 0 (bit-and y 0x1))
                0
                matrix-a))))
                
(defn -genrand-nth []
	 "Returns the next number from the seuquence (factored from the reference implementation)"
   (dosync
    (if 
    	;; Do we need to regenerate the state?
      (>= @mti n)
        (do
          (if (= @mti (inc n))
            (sgenrand 4357))
          ;; Build a new vector of integers from the old one (and itself!)
          (ref-set mt  (loop [k (int 0) v []]
                          (if (< k n)
                             (let [v-1 (cond 
                                         (< k (- n m)) (conj v (-genrand-inner (nth @mt k) (nth @mt (inc k)) (nth @mt (+ k m))))
                                         (< k (dec n)) (conj v (-genrand-inner (nth @mt k) (nth @mt (inc k)) (nth v (+ k (- m n)))))
                                         (= k (dec n)) (conj v (-genrand-inner (nth @mt (dec n)) (nth v 0) (nth v (dec m)))))]
                             (recur (inc k) v-1))
                            v)))
           (ref-set mti 0)))
      ;; Now return the element and advance the counter
      (let [retval (int (nth @mt @mti))]
        (ref-set mti (inc @mti))
        retval)))
   
(defn genrand []
	"Generate the next number in the sequence, if this number is the nth in the sequence then
	 a new state table will be generated in this call"	
  (let [y (-genrand-nth)
        y1 (bit-xor y (tempering-shift-u y))
        y2 (bit-xor y1 (bit-and (tempering-shift-s y1) tempering-mask-b))
        y3 (bit-xor y2 (bit-and (tempering-shift-t y2) tempering-mask-c))
        y4 (bit-xor y3 (tempering-shift-l y3))]
     y4))
      
