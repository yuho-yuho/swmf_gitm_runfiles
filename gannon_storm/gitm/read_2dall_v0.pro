PRO read_2dall

  SAT_KEY='2dall'

  DATA_DIR = './'

  SAVE_DIR ='./'
  SAVE_DIR = SAVE_DIR + SAT_KEY +'/'

  SAVE_DIR_EXIST=FILE_TEST(SAVE_DIR,/DIRECTORY)                               
                                                                               
  IF (SAVE_DIR_EXIST gt 0) THEN BEGIN                                          
     print, 'SAVE DIR EXISTS !'                                                
                                                                               
  ENDIF ELSE BEGIN                                                             
                                                                               
     FILE_MKDIR, SAVE_DIR                                                      
     print, 'SAVE DIR HAS BEEN CREATED !'                                      
                                                                               
  ENDELSE

  FLIST=FINDFILE(DATA_DIR+'2DUSR*.bin')
  NFILES=N_ELEMENTS(FLIST)                                                     
                                                                               
  print, nfiles

  ;nfiles=1

  FOR IFILE=0,nfiles-1 DO BEGIN                                                
                                                                               
                                                                               
     FILENAME = FLIST(IFILE)                                                   
     PRINT, IFILE;,FILENAME      
                                                                               
     data = 0.0                                                                
                                                                               
     read_thermosphere_file, filename, nvars, nalts, nlats, nlons, $           
                             vars, data, rb, cb, bl_cnt, itime1                
                                                                               
                                                                               
     ;alt = reform(data(2,*,*,*))                      
     lat = reform(data(1,*,*,*))                                 
     lon = reform(data(0,*,*,*))                        
     
     ;help,alt,lat,lon

     ;for i=0,nvars-1 do print, tostr(i)+'. '+vars(i) 
  
     Potential = reform(data(4,*,*,*))
     Tec = reform(data(5,*,*,*))
     intJH = reform(data(6,*,*,*))
     Eflux = reform(data(3,*,*,*))   
     jp = reform(data(14,*,*,*))
     jwd = reform(data(15,*,*,*))
     jef = reform(data(16,*,*,*))   

     SAVE_STR=strmid(FILENAME,0,20)
     print, save_str
     SAVE_FN=SAVE_DIR+save_str+'.sav'

     SAVE, itime1, lon, lat,  Potential, Eflux, intJH, $
           Tec,  jp, jwd, jef, FILENAME=SAVE_FN

  ENDFOR

END
