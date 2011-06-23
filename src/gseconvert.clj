(ns gseconvert
  (:require
   [clojure.contrib.io :as io]))

(defn ensure-unzip [file]
  ;;TODO: implement
  file)

(defn read-soft-section [file section]
  (take-while #(not (.startsWith % (str "!" section "_end")))
              (rest
               (drop-while
                #(not (.startsWith % (str "!" section "_begin")))
                (io/read-lines
                 (ensure-unzip file))))))

;;Functions to read SOFT platform matrix and map probes to Entrez Gene IDs
(defn default-mapper [m k]
  (zipmap (:ID m)
          (map (fn [g]
                (if-not (empty? g)
                  (map #(Integer/parseInt (.trim %))
                       (.split g "///"))))
               (m k))))

(def make-mapper nil)
(defmulti make-mapper
  #(first
   (filter identity
           (map #{:ENTREZ_GENE_ID}
                (keys %)))))

(defmethod make-mapper :ENTREZ_GENE_ID [m]
  (default-mapper m :ENTREZ_GENE_ID))

(defn probes-to-entrez-ids* [file]
  "Read platform specification file (e.g. 'GPL96.txt'), and return
a map of probes to Entrez Gene IDs, if possible, otherwise nil."
  (make-mapper
   (let [lines (read-soft-section file "platform_table")
         ks (.split (first lines) "\t")]
     (zipmap (map keyword ks)
             (apply map vector
                    (map #(let [fields (.split % "\t")]
                            (concat fields (repeat (- (count ks) (count fields)) nil)))
                         (take-while #(not (.startsWith % "!platform_table_end"))
                                     (rest lines))))))))
(defn probes-to-entrez-ids [gpl]
  ;;TODO: implement
  (probes-to-entrez-ids*
   "/home/gilesc/Desktop/GPL96.txt"))


(defn mean [coll]
  (/ (apply + coll)
     (count coll)))

(defn scale [expression]
  )

(defn scale-and-validate [expression-vec]
  "Various postprocessing steps:
   - Floor 'extreme values': make the top 0.1% of values equal to the minimum of those top 0.1%.

   - Scale each experiment to 0-10000.
  
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
      (let [max-value (nth (reverse (sort expression-vec))
                           (math/ceil (/ (count expression-vec) 1000)))
            min-value (apply min expression-vec)]
        ;;;;;THIS IS WHERE I STOPPED FIX THIS SCALING
        ))))

(defn read-expression [file gpl]
  (let [mapper (probes-to-entrez-ids gpl)
        lines (read-soft-section file "series_matrix_table")
        gsms (rest (.split (first lines) "\t"))
        mtx ;;Take signature with highest average expression
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
    {:row-names gsms :col-names genes :expression mtx-t}))

(defn into-global-expression-vector [expression])
