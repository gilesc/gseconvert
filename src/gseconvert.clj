(ns gseconvert
  (:require
   [clojure.contrib.math :as math]
   [clojure.contrib.io :as io]
   [clojure.contrib.string :as string]
   [clojure.contrib.shell-out :as shell]))

(defn ensure-unzip [file]
  (if (.endsWith (str file) ".gz")
    (java.util.zip.GZIPInputStream.
     (java.io.FileInputStream. file))
    file))

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

(defn get-gpl-file [gpl]
  (let [base (str "data/GPL/" gpl ".annot.gz")]
    (first (filter #(.exists %)
                   (map #(java.io.File. %) [(str base ".1") base])))))

(defn probes-to-entrez-ids [gpl]
  "Read platform specification file (e.g. 'GPL96.txt'), and return
a map of probes to Entrez Gene IDs, if possible, otherwise nil."
  (if-let [file (get-gpl-file gpl)]
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
    {:row-names gsms :col-names genes :expression mtx-t}))

(require '[clojure.contrib.pprint :as pprint])

(defn pprint-head [e]
  (pprint/pprint (map #(take 5 %) (take 5 (:expression e)))))

(defn into-global-expression-vector [genes expression]
  (for [row (:expression expression)
        :let [m (zipmap (:col-names expression) row)]]
    (map #(m % "NA") genes)))

(defn shell-exec [cmd]
  (shell/sh "sh" "-c" cmd))

(defn get-tax-id [species]
  (Integer/parseInt
   (.trim
    (shell-exec
     (format "grep -P \"\t%s\t\" data/taxonomy.dat | cut -f1 | head -1" species)))))

(defn get-genes [tax-id]
  (map #(Integer/parseInt %)
       (string/split #"\n"
        (shell-exec
         (format "grep -P \"^%s\t\" data/gene_info | cut -f2" tax-id)))))

(defn -main [species]
  (prn species))
