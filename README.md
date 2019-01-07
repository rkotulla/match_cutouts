# match_cutouts

This is a small repo to keep track of a tool to match images based on header data. it uses a reference frame (either fuill frame or just header) and then matches a second FITS image to this reference fraem. 

Note that this only works if the reference frame does not include any distortions. If it does you need to run it though swarp first to account for distortions, and then you can use this tool to match the other frame.
