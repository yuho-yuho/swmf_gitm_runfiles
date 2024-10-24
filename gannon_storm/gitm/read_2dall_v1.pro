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
  
     dbt_e = reform(data(26,*,*,*))
     dbt_n = reform(data(27,*,*,*))
     dbt_z = reform(data(28,*,*,*))
     dbwd_e = reform(data(29,*,*,*))   
     dbwd_n = reform(data(30,*,*,*))
     dbwd_z = reform(data(31,*,*,*))
     dbef_e = reform(data(32,*,*,*))
     dbef_n = reform(data(33,*,*,*))
     dbef_z = reform(data(34,*,*,*))

     SAVE_STR=strmid(FILENAME,0,20)
     print, save_str
     SAVE_FN=SAVE_DIR+save_str+'.sav'

     SAVE, itime1, lon, lat,  dbt_e, dbt_n, dbt_z, dbwd_e,  $
           dbwd_n, dbwd_z, dbef_e, dbef_n, dbef_z, FILENAME=SAVE_FN

  ENDFOR

END
