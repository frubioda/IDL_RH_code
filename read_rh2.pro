; *******************************************************************

pro read_rh2,dir,spect,atoms,atmosf,geom

@input.common
@geometry.common
@files.common
@atmos.common
@spectrum.common
@opacity.common

  cd,dir,current=calling_dir
  
  print, "-- Reading inputData ...."
  IF (readInput('input.out')) THEN print, "     Success" ELSE print, "     Failed"
  
  print, "-- Reading geometry ...."
  IF (readgeometry('geometry.out')) THEN print, "     Success" ELSE print, "     Failed"
  
geom=geometry

  print, "-- Reading atmosphere ...."
  metalFile  = "metals.out"  &  moleculeFile = "molecules.out"
  IF (readatmos('atmos.out')) THEN print, "     Success" ELSE print, "     Failed"
  
atmosf=atmos

  print, "-- Reading spectrum ...."
  IF (readspectrum('spectrum.out')) THEN print, "     Success" ELSE print, "     Failed"
  
  spect=spectrum

  ;; print, "-- Reading atomic models ...."
  ;; atom_files = file_search("atom.*.out", COUNT=Nfile)
  ;; print, "Found ", Nfile, " atoms: "
  ;; IF (Nfile GT 0) THEN print, "  atom_files: ", atom_files
  ;; atom=readatom(atom_files[0])
  ;; natom=n_elements(atom_files)
  ;; atoms=replicate(atom,natom)
  ;; for i=0,natom-1 do atoms[i]=readatom(atom_files[i])
  
  cd,calling_dir

end
