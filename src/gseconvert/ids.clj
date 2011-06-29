(ns gseconvert.ids
  "Module responsible for converting among identifiers, such as GSM, GPL, and GSE,
or taxonomy string, taxonomy id, and Entrez Gene IDs."
  (:require
   [clojure.contrib.sql :as sql]
   [clojure.contrib.shell-out :as shell]
   [clojure.contrib.string :as string]))

(let [db {:classname "org.sqlite.JDBC"
          :subprotocol "sqlite"
          :subname "data/GEOmetadb.sqlite"}]
  (defn- query-db [qry]
    (sql/with-connection db
      (sql/with-query-results rs
        [qry]
        (doall rs)))))

(defn shell-exec [cmd]
  (shell/sh "sh" "-c" cmd))


(def convert-id nil)
(defmulti convert-id
  "Use like: (convert-id [:gpl :gse] \"GPL96\"), where the first argument is a vector
[:from :to]"

  (fn [from-to _] from-to))

(defmethod convert-id [:gpl :gse] [_ gpl]
  {:pre [(.startsWith gpl "GPL")]}
  (map :gse
   (query-db
    (format "SELECT DISTINCT(gse) FROM gse_gpl WHERE gpl='%s'" gpl))))

(defmethod convert-id [:species :gpl] [_ species]
  "Returns GPLs associated with a given species string.
Additionally guarantees that GPLs are sorted from most to least common."
  {:pre [(string? species)]}
  (map :gpl
   (reverse
    (sort-by :count
             (query-db
              (format "SELECT gpl,count(gpl) as count FROM gsm WHERE organism_ch1='%s' GROUP BY gpl" species))))))

(defmethod convert-id [:species :tax-id] [_ species]
  (map #(Integer/parseInt %)
       (string/split #"\n"
        (shell-exec
         (format "grep -P \"\t%s\t\" data/taxonomy.dat | cut -f1" species)))))

(defmethod convert-id [:tax-id :gene] [_ tax-id]
  (map #(Integer/parseInt %)
       (string/split #"\n"
        (shell-exec
         (format "grep -P \"^%s\t\" data/gene_info | cut -f2" tax-id)))))

(defmethod convert-id [:species :gene] [_ species]
  (convert-id [:tax-id :gene]
              (first
               (convert-id [:species :tax-id] species))))

(defmethod convert-id [:gpl :species] [_ gpl]
  (map :organism
   (query-db
    (format "SELECT organism FROM gpl WHERE gpl='%s'" gpl))))

(defmethod convert-id [:gpl :gene] [_ gpl]
  (convert-id [:species :gene]
   (first (convert-id [:gpl :species] gpl))))
