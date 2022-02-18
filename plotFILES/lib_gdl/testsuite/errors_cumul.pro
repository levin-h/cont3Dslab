;
; Alain C., 28 Feb. 2017
;
; Feb. 2018 : Change in naming convention : 
; ERRORS_CUMUL, ERRORS_ADD, ERRORS_RESET
;
; -----------------------------------------------
;
; Purpose : managing Errors Accumulation
; This procedure adds running "new_errors" into "cumul_errors"
;
; It would surprising if "new_errors" is undefined
; It is *not* surprising that "cumul_errors" may be undefined
; 
; -----------------------------------------------
; In order to test fake errors accumulations in the tests, 
; you can use : ERRORS_CUMUL, /debug 
;
; BUT we will have to change it for a better ON/Off triggering
;
; -----------------------------------------------
;
pro ERRORS_CUMUL, cumul_errors, new_errors, help=help, $
                  debug=debug, verbose=verbose, test=test
;
; this case /debug should be call BEFORE or OUTSIDE
;
if KEYWORD_SET(debug) then begin
    DEFSYSV, '!cumul', 1
    MESSAGE, /continue, 'Fake error has be added now'
    return
endif
;
if KEYWORD_SET(help) OR (N_PARAMS() NE 2) then begin
    print, 'Usage : pro ERRORS_CUMUL, cumul_errors, new_errors, help=help, $'
    print, '                          debug=debug, verbose=verbose, test=test'
    return
endif
;
if KEYWORD_SET(verbose) then begin
    print, 'Value of >>cumul_errors<< : ', cumul_errors
    print, 'Value of >>new_errors<<   : ', new_errors
endif
;
if (SIZE(cumul_errors, /type) GT 0) then begin
    if (SIZE(new_errors, /type) GT 0) then begin
        cumul_errors=cumul_errors+new_errors
    endif
    ;; nothing to add if "new_errors" not defined
endif else begin
    if (SIZE(new_errors, /type) GT 0) then begin
        cumul_errors=new_errors
    endif
endelse
;
if KEYWORD_SET(verbose) then begin
    print, 'NEW Value of >>cumul_errors<< : ', cumul_errors
endif
;
; debug mode
DEFSYSV, '!cumul', exist=exist
if (exist) then begin
   if (!cumul GT 0) then begin
      if (SIZE(cumul_errors,/type) GT 0) then cumul_errors++ else cumul_errors=1
   endif
endif
;
if KEYWORD_SET(verbose) then begin
    print, 'NEW Value of >>cumul_errors<< : ', cumul_errors
endif
;
if KEYWORD_SET(test) then STOP
;
end
