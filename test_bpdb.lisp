;; #!/usr/bin/env sbcl --script

(require 'asdf)

;; --- 1. Utilities ---

(defun starts-with-p (prefix string)
  (let ((plen (length prefix)))
    (and (>= (length string) plen)
         (string= prefix string :end2 plen))))

(defun trim-r (string)
  (string-right-trim '(#\Space #\Tab #\Newline #\Return) string))

;; --- 2. Normalizer ---

(defun process-atom-line (line)
  (let* ((len (length line))
         (padded (if (< len 80)
                     (concatenate 'string line (make-string (- 80 len) :initial-element #\Space))
                     line)))
    (when (starts-with-p "HETATM" padded)
      (setf (subseq padded 0 6) "ATOM  "))
    (when (>= (length padded) 11) (setf (subseq padded 6 11) "     "))
    (when (>= (length padded) 17) (setf (char padded 16) #\Space))
    (when (and (>= (length padded) 60) (char= (char padded 56) #\1))
      (setf (subseq padded 57 60) ".00"))
    (let ((chg (subseq padded 78 80)))
      (when (or (string= chg "0 ") (string= chg " 0") (string= chg "+0") 
                (string= chg "-0") (string= chg "0+") (string= chg "0-") (string= chg "00"))
        (setf (subseq padded 78 80) "  ")))
    (trim-r padded)))

(defun supported-p (line)
  (let ((prefixes '("HEADER" "TITLE " "EXPDTA" 
                    "SSBOND" "CISPEP" "HELIX " "SHEET " "SITE  "
                    "ATOM  " "HETATM" "TER   ")))
    (loop for p in prefixes thereis (starts-with-p p line))))

(defun normalize-stream (in-stream)
  (let ((res-blocks nil)
        (current-block nil)
        (current-key ""))
    (loop for line = (read-line in-stream nil nil)
          while line do
      (let ((norm (trim-r line)))
        (when (supported-p norm)
          (cond
            ((or (starts-with-p "ATOM  " norm) (starts-with-p "HETATM" norm))
             (let* ((processed (process-atom-line norm))
                    (padded (if (< (length processed) 80)
                                (concatenate 'string processed (make-string (- 80 (length processed)) :initial-element #\Space))
                                processed))
                    (name (subseq padded 12 16))
                    (elem (trim-r (string-left-trim '(#\Space) (subseq padded 76 78)))))
               (unless (or (find #\H name) (string= elem "H") (string= elem "D"))
                 (let ((key (subseq padded 21 27)))
                   (if (string= key current-key)
                       (push processed current-block)
                       (progn
                         (when current-block
                           (push (stable-sort current-block #'string< 
                                   :key (lambda (x) (subseq (concatenate 'string x "                    ") 12 16))) res-blocks))
                         (setf current-key key)
                         (setf current-block (list processed))))))))
            ((starts-with-p "TER   " norm)
             (when current-block
               (push (stable-sort current-block #'string< 
                       :key (lambda (x) (subseq (concatenate 'string x "                    ") 12 16))) res-blocks)
               (setf current-block nil)
               (setf current-key ""))
             (push (list (process-atom-line norm)) res-blocks))
            ((starts-with-p "SHEET " norm)
             (when current-block
               (push (stable-sort current-block #'string< 
                       :key (lambda (x) (subseq (concatenate 'string x "                    ") 12 16))) res-blocks)
               (setf current-block nil)
               (setf current-key ""))
             (if (> (length norm) 38)
                 (push (list (trim-r (subseq norm 0 38))) res-blocks)
                 (push (list norm) res-blocks)))
            (t
             (when current-block
               (push (stable-sort current-block #'string< 
                       :key (lambda (x) (subseq (concatenate 'string x "                    ") 12 16))) res-blocks)
               (setf current-block nil)
               (setf current-key ""))
             (push (list norm) res-blocks))))))
    (when current-block
      (push (stable-sort current-block #'string< 
              :key (lambda (x) (subseq (concatenate 'string x "                    ") 12 16))) res-blocks))
    (loop for block in (nreverse res-blocks)
          append block)))

(defun normalize-file (path)
  (with-open-file (in path :direction :input :external-format :utf-8)
    (normalize-stream in)))

(defun normalize-string (str)
  (with-input-from-string (in str)
    (normalize-stream in)))

(defun compare-lists (l1 l2)
  (loop for x in l1
        for y in l2
        for i from 1
        unless (string= x y)
        do (return (format nil "Mismatch at line ~D:~%Orig: '~A'~%Dec : '~A'" i x y))
        finally (if (= (length l1) (length l2))
                    (return nil)
                    (return (format nil "Length mismatch: Original ~D lines, Decoded ~D lines." (length l1) (length l2))))))

;; --- 3. Test Runner ---

(defun run-test-on-file (pdb-file bpdb-bin)
  (let* ((orig-lines (normalize-file pdb-file))
         (cmd (format nil "~A encode ~A - | ~A test - -" bpdb-bin (namestring pdb-file) bpdb-bin))
         (decoded-str (uiop:run-program (list "sh" "-c" cmd) :output :string :ignore-error-status t))
         (dec-lines (normalize-string decoded-str)))
    (let ((diff (compare-lists orig-lines dec-lines)))
      (if diff
          (values nil diff)
          (values t nil)))))

;; --- 4. File Caching Utility ---

(defun get-pdb-files (db-dir wildcard)
  (let ((cache-file "pdb_files_cache.txt"))
    (if (probe-file cache-file)
        (progn
          (format t "Reading file list from cache: ~A~%" cache-file)
          (with-open-file (in cache-file :direction :input)
            (loop for line = (read-line in nil nil)
                  while line
                  collect (parse-namestring line))))
        (progn
          (format t "Scanning ~A for .pdb files (this may take a while)...~%" db-dir)
          (let ((files (directory wildcard)))
            (format t "Caching ~D files to ~A~%" (length files) cache-file)
            (with-open-file (out cache-file :direction :output :if-exists :supersede)
              (dolist (f files)
                (format out "~A~%" (namestring f))))
            files)))))

;; --- 5. CLI Execution ---

(defun main ()
  (let ((args (cdr sb-ext:*posix-argv*)))
    (when (< (length args) 1)
      (format t "Usage: ./test-suite.lisp <path-to-pdb-dir> [num-samples]~%")
      (uiop:quit 1))
    
    (let* ((db-dir (first args))
           (num-samples (if (> (length args) 1) (parse-integer (second args)) 1000))
           (bpdb-bin "./bpdb")
           (wildcard (make-pathname :directory (append (pathname-directory (truename db-dir)) 
                                                       (list :wild-inferiors))
                                    :name :wild :type "pdb")))
      
      (setf *random-state* (sb-ext:seed-random-state 1238))
      
      (with-open-file (log-stream "test_results.log" :direction :output 
                                  :if-exists :supersede :if-does-not-exist :create)
        (let ((all-files (get-pdb-files db-dir wildcard)))
          (unless all-files
            (format t "No .pdb files found.~%")
            (uiop:quit 1))
          
          (format t "Found ~D files. Selecting ~D with seed 42.~%" (length all-files) num-samples)
          (format log-stream "Found ~D files. Selecting ~D with seed 42.~%" (length all-files) num-samples)
          
          (let ((file-vec (coerce all-files 'vector)))
            (loop for i from 0 below (length file-vec)
                  for j = (+ i (random (- (length file-vec) i)))
                  do (rotatef (aref file-vec i) (aref file-vec j)))
            (setf all-files (coerce file-vec 'list)))
            
          (let ((samples (subseq all-files 0 (min num-samples (length all-files))))
                (passed 0)
                (failed 0))
            
            (loop for f in samples
                  for i from 1 do
              (multiple-value-bind (ok diff-out) (run-test-on-file f bpdb-bin)
                (if ok
                    (progn
                      (incf passed)
                      (format t "[~D/~D] [PASS] ~A~%" i (length samples) (pathname-name f))
                      (format log-stream "[~D/~D] [PASS] ~A~%" i (length samples) (pathname-name f)))
                    (progn
                      (incf failed)
                      (format t "[~D/~D] [FAIL] ~A (see log)~%" i (length samples) (pathname-name f))
                      (format log-stream "[~D/~D] [FAIL] ~A~%~A~%" i (length samples) (pathname-name f) diff-out)))))
            
            (let ((summary (format nil "~%=== SUMMARY ===~%Passed: ~D~%Failed: ~D~%" passed failed)))
              (format t "~A" summary)
              (format log-stream "~A" summary)
              (if (= failed 0)
                  (progn
                    (format t "🎉 All structures match perfectly!~%")
                    (format log-stream "🎉 All structures match perfectly!~%"))
                  (progn
                    (format t "⚠️ Check test_results.log for diffs.~%")
                    (format log-stream "⚠️ Check test_results.log for diffs.~%"))))))))))

(main)
