(ns gseconvert.files
  "Given GPL and GSE identifiers, resolve and return the correct file(s)."
  (:require
   [clojure.contrib.io :as io]))

(defn ensure-unzip [file]
  (if (.endsWith (str file) ".gz")
    (java.io.InputStreamReader.
     (java.util.zip.GZIPInputStream.
      (java.io.FileInputStream. file)))
    file))

(defn existing-file [& paths]
  "Return the first file from a list of (String) paths that actually exists."
  (first
   (filter #(.exists %)
           (map #(java.io.File. %)
                paths))))

(defn read-gse [gse gpl]
  (let [base (str "data/GSE/" gse)
        suffix "_series_matrix.txt.gz"]
    (ensure-unzip
     (existing-file (str base "-" gpl suffix)
                    (str base suffix)))))

(defn read-gpl [gpl]
  (let [base (str "data/GPL/" gpl ".annot.gz")]
    (ensure-unzip
     (existing-file base (str base ".1")))))

(defn read-soft-section [file section]
  (take-while #(not (.startsWith % (str "!" section "_end")))
              (rest
               (drop-while
                #(not (.startsWith % (str "!" section "_begin")))
                (io/read-lines file)))))

