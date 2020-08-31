function ProfileScript(script_name, profile_dir)
  ## ProfileScript - Profile a script.
  ## @param script_name The script name for profiling.

  if (!exist("profile_dir", "var"))
    profile_dir = cstrcat("./", script_name);
  endif

  ## Start profiling for the algorithm execution.
  profile clear;
  profile on;

  eval(script_name);

  ## Stop the profiling.
  profile off;
  ## Get profile data.
  prof_data = profile("info");
  ## Show the profile data in a flat view.
  profshow(prof_data, length(prof_data.FunctionTable));
  ## Export the profile data to HTML files.
  profexport(profile_dir);
endfunction
