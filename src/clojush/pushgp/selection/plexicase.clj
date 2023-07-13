(ns clojush.pushgp.selection.plexicase
  (:require (bigml.sampling [simple :as simple])))

(def arithmetic-mean
  "Calculates the arithmetic mean of numeric inputs."
  (fn
    ([] [0.0 0])
    ([[^double s ^long c :as acc] e]
     (if (nil? e)
       acc
       (let [e (double e)]
         [(+ s e) (inc c)])))
    ([[s c]]
     (when-not (zero? c)
       (/ s c)))))

(defn mean
  "Arithmetic mean."
  [arr]
  (transduce identity arithmetic-mean arr))

(defn calculate-elitist-on-case
  "Calculates the hj value for the individual's error on case and the population
   min error on case."
  [ind-error pop-min-error]
  (if (= ind-error pop-min-error) 1 0))

(defn calculate-hj-vector
  "Calculates the hj values for this individual."
  [individual min-error-per-case option]
  (let [ind-elitist-vector (mapv calculate-elitist-on-case
                                (:errors individual)
                                min-error-per-case)]
    (case option
      1 (let [;; E_yi is the number of cases on which this individual is elite
              E_yi (apply + ind-elitist-vector)]
          (replace {1 E_yi} ind-elitist-vector))
      2 ind-elitist-vector
      (throw (Exception. (str "Unrecognized plexicase option: " option))))))

(defn unnormalized-probability-distribution
  "hj function from plexicase paper across whole population
   Calculate probability distribution for individuals across all test cases

   Options: 
   individuals: a vector of individuals 
   num-cases: the number of test cases

   Returns: 
   A list of lists of unnormalized probabilities for each individual"
  [{:keys [individuals num-cases option]}]
  (let [min-error-per-case (mapv (fn [c] (reduce min (mapv #(nth (:errors %) c)
                                                         individuals)))
                                (range num-cases))]
    (mapv #(calculate-hj-vector % min-error-per-case option)
         individuals)))
     
(defn normalize-probability-distribution
  "Calculates the Pj(yi) values for every case and every individual."
  [unnormalized-prob-distribution]
  (let [row-sums (mapv #(reduce + %) (apply mapv list unnormalized-prob-distribution))]
    (mapv (fn [hj-vector]
           (mapv (fn nested-normalize-prob-dist-anon [hj row-sum] (float (/ hj row-sum)))
                hj-vector
                row-sums))
         unnormalized-prob-distribution)))

(defn probability-distribution
  "Calculates P(yi) = the final probability for each individual."
  [normalized-prob-dist]
  (mapv mean normalized-prob-dist))

(defn plexicase-select-all-parents
  "Selects all parents for the generation. This is done here instead of one
   parent at a time, because simple/sample is much faster that way."
  [individuals probability-distribution number-parents]
  (vec(take number-parents
        (map #(nth individuals %)
             (simple/sample (range (count individuals))
                            :weigh (fn [indi] (nth probability-distribution indi))
                            :replace true)))))

(defn make-plexicase-selection
  "Calculates the selection probabilities for every individual in the population.
   Params:
     state - a map that contains at least :individuals as a key. Can optionally
             take option, which is a number indicating which version of plexicase
             to use:
               option 1: original plexicase
               option 2: uses 1 instead of E(yi) in the calculation of hj"
  [num-parents {:keys [individuals option num-errors] :or {option 1}}]
  (let [num-cases num-errors
        unnormalized-prob-dist (unnormalized-probability-distribution {:individuals individuals :num-cases num-cases :option option})
        normalized-prob-dist (normalize-probability-distribution unnormalized-prob-dist)
        ind-probabilities (probability-distribution normalized-prob-dist)]
    (plexicase-select-all-parents individuals ind-probabilities num-parents)))

(defn plexicase-select-parent-using-index
  "Uses plexicase to select a parent. Uses the parents already selected by
   plexicase-select-all-parents"
  [plexicase-parents index]
  (nth plexicase-parents index))
