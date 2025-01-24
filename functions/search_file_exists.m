function file_search_indx = search_file_exists(dir_path, file_name)

folder_files = dir(dir_path);

file_names = {folder_files.name};

file_search_indx = strcmp(file_names, file_name);

end