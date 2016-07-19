def py_peak_calling(bedgraph, threshold, min_length, inter_peak_distance, max_length=10000,
                   generate_ID=True, output_name = None):
    """
    Written by Pete Skene (peteskene@gmail.com). Free for academic use.
    
    - need to install a more up-to-date varsion of bedtools before invoking Jupyter
      type: module load bedtools/2.21.0
    - (1) filters bedgraph based on threshold; (2) merges adjacent basepairs that are over threshold;
      (3) retains peaks that satisfy min/max length criteria; (4) merges any peaks that are closer
      than the inter-peak distance cutoff
    - max length is typically defaulted to be very large
    - outputs a bed file (default col4 is the sum of the bedgraph scores; sorted by chrom;start;stop)
    - generate ID: will auto generate a integer list as a ID number (1... number of peaks). This will 
    be reported as column 4 and the bedgraph scores will be shifted to column 5 as per standard bed format
    - note the peak score for merged peak is the *just* the sum of the two individual peaks not the 
    total score in the merged region (i.e. there could be some sub-threshold scores in the intervening 
    space that won't be included)
    -assumes bedgraph in standard format <chr> <start> <stop> <score>
    -output_name = option for user defined name (type with '...'), otherwise will generate name bedgraph_peaks.bed
    """
    
    import pybedtools
    import glob
    from pybedtools import BedTool
    import pandas as pd
    
    #generate name for output
    bedgraph_name = glob.glob(bedgraph)
    
    if output_name != None:
        filename = output_name
        
    elif output_name == None:
        filename = bedgraph_name[0].replace('.bg', '_peaks.bed')
        
    print 'input bedgraph file: ' + bedgraph_name[0]
    print 'output filename: ' + filename
    
    #import data as BedTool
    data = BedTool(bedgraph) 
    
    #retains intervals above threshold
    above_thresh = data.filter(lambda b: float(b.name) >= threshold) 
    
    #merge adjacent above threshold regions and sum bedgraph scores (assumes bedgraph score in col 4)
    #by increasing d value can allow for 
    merge_regions= above_thresh.merge(d=0, c=4, o='sum' )
    
    #filter based on length criteria
    peaks = BedTool(merge_regions.filter(lambda x: len(x) >= min_length and len(x) <= max_length))
    
    #merge the bonafide peaks if they they are shorter than the inter peak distance and sum scores and sort
    merge_peaks = peaks.merge(d=inter_peak_distance, c= 4, o='sum').sort()
    
    print 'number of peaks found: ' + str(merge_peaks.count())
    
    if not generate_ID:
        print 'saving sorted peak bed file with no ID'
        
        merge_peaks.saveas(filename)
        
    if generate_ID:
        print 'saving sorted peak bed file with ID names'
        
        #change to pandas dataframe
        DF_peaks = merge_peaks.to_dataframe()
        
        #insert new column with id: 1.... # of peaks
        DF_peaks.insert(3, 'id', ['id' + str(item) for item in range(1, (len(DF_peaks)+1))])
        
        ['id' + str(item)  for item in range(1, 5)]
        #save output
        DF_peaks.to_csv(filename, sep = '\t', header = False, index = False)
        
    return 'Finished'
    
