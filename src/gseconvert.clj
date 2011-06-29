(ns gseconvert
  (:use
   gseconvert.files
   gseconvert.ids)
  (:require
   [clojure.contrib.math :as math]
   [clojure.contrib.io :as io]
   [clojure.contrib.string :as string]))

;;Functions to read SOFT platform matrix and map probes to Entrez Gene IDs
(defn default-mapper [m k]
  (zipmap (:ID m)
          (map (fn [g]
                (if-not (empty? g)
                  (map #(Integer/parseInt (.trim %))
                       (.split g "///"))))
               (m k))))

(defn make-mapper [m]
  (let [k (first
           (filter identity
                   (map #{:ENTREZ_GENE_ID (keyword "Gene ID")}
                        (keys m))))]
    (default-mapper m k)))


(defn probes-to-entrez-ids [gpl]
  "Read platform specification file (e.g. 'GPL96.txt'), and return
a map of probes to Entrez Gene IDs, if possible, otherwise nil."
  (if-let [file (read-gpl gpl)]
    (make-mapper
     (let [lines (read-soft-section file "platform_table")
           ks (.split (first lines) "\t")]
       (zipmap (map keyword ks)
               (apply map vector
                      (map #(let [fields (.split % "\t")]
                              (concat fields (repeat (- (count ks) (count fields)) nil)))
                           (take-while #(not (.startsWith % "!platform_table_end"))
                                       (rest lines)))))))))

(defn mean [coll]
  (/ (apply + coll)
     (count coll)))

(defn median [coll]
  ((vec (sort coll))
   (math/round (/ (count coll) 2))))


(defn scale [expression-vec]
  "Scale each experiment to 0-10000.
   Floor 'extreme values': make the top 0.1% of values equal to the minimum of those top 0.1%. "
  (let [max-value (nth (reverse (sort expression-vec))
                           (math/ceil (/ (count expression-vec) 1000)))
        min-value (apply min expression-vec)]
        (map #(* 10000
                 (/ (- % min-value)
                    (max 1 (- max-value min-value))))
             expression-vec)))

(defn scale-and-validate [expression-vec]
  "Various postprocessing steps:
   - Scaling (see above).
   - Detect and fix experiments with logratio instead of raw expression values. (TODO)

   Also, validate each experiment according to criteria:
      1. Mean and median >= 0
      2. Mean-median ratio >= 1.2
      3. <= 1% negative values

   Returns nil for experiments not meeting all criteria."
  ;;TODO: detect log experiments
  (let [e-mean (mean expression-vec)
        e-median (median expression-vec)]
    (if (and (>= e-mean 0)
             (>= e-median 0)
             (>= (/ e-mean e-median) 1.2)
             (/ (count (filter neg? expression-vec)) (count expression-vec)))
      (scale expression-vec))))

(require '[clojure.contrib.pprint :as pprint])
(defn pprint-head [e] ;;TODO: rm
  (pprint/pprint (map #(take 5 %) (take 5 (:expression e)))))

(defn read-expression [file gpl]
  (let [mapper (probes-to-entrez-ids gpl)
        lines (read-soft-section file "series_matrix_table")
        gsms (rest (.split (first lines) "\t"))
        mtx ;;Takes signature with highest average expression
        (into {}
         (for [[gene signatures]
               (group-by first
                         (apply concat
                                (for [line (rest lines)
                                      :let [fields (.split line "\t")
                                            expression
                                            (map #(Float/parseFloat %)
                                                 (rest fields))]]
                                  (for [gene (mapper (.replaceAll (first fields) "\"" ""))]
                                    [gene expression]))))]
           [gene (last (sort-by #(- (mean %)) (map second signatures)))]))
        genes (keys mtx)
        mtx-t (filter identity
                      (map scale-and-validate
                           (apply map vector (vals mtx))))]
    (if-not (empty? (first mtx-t))
      {:row-names gsms :col-names genes :expression mtx-t})))

(require '[clojure.contrib.pprint :as pprint])


(defn to-global-expression-vectors [genes expression-mtx]
  (prn (count (:expression expression-mtx)))
  "Creates one column for every gene in the species and returns a Seq of ordered vectors."
  (map #(cons %1 %2)
       (:row-names expression-mtx)
       (for [row (:expression expression-mtx)
             :let [m (zipmap (:col-names expression-mtx) row)]]
         (map #(m % "NA") genes))))

(defn expression-matrix [gpl]
  (let [genes (sort (convert-id [:gpl :gene] gpl))]
    (cons
     (concat ["-"] genes)
     (for [gse (take 5 (convert-id [:gpl :gse] gpl))
           row (to-global-expression-vectors genes
                                                 (read-expression (read-gse gse gpl) gpl))]
       row))))

(defn write-tsv [file rows]
  (io/write-lines file
   (map #(string/join "\t" %)
        rows)))

(defn species-expression-matrix [species]
  (apply concat
         (map expression-matrix
              (take 3
                    (convert-id [:species :gpl] species)))))

(defn -main [species]
  (write-tsv (format "data/gse-%s.mtx"
                     (.replaceAll (.toLowerCase species) " " "_"))
             (species-expression-matrix species)))
