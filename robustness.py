#################################################
#   robustness.py                               #
#   author: Byron C. Wallace                    #
#                                               #
#   This module contains code for our 'How      #
#   sensitive are meta-analyses to missed       #
#   articles' paper.                            #
#################################################

import pdb
import random
import math
import sys

try:
    import rpy2
    from rpy2 import robjects as ro 
    ro.r("library(metafor)")
except:
    print "Whoops -- there was a problem import rpy2 -- is this module installed?"

####################################
# Objects to map the database to   #
####################################
BINARY, CONTINUOUS = range(2)

### path to .tsv file containing binary data.
binary_data_path = "binary_database_with_method.tsv"

# maps indices (magic numbers) to strings.
methods_d = {1: 'peto OR',
             2: 'OR Fixed (MH)',
             3: 'OR Random (based on MH Fixed)',
             4: 'RR Fixed (MH)',
             5: 'RR Random (based on MH Fixed)',
             6: 'RD Fixed (MH)',
             7: 'RD Random (based on MH Fixed)',
             8: 'WMD Fixed',
             9: 'WMD Random',
             10: 'SMD Fixed',
             11: 'SMD Random'}

class Review:
    ''' Reviews (potentially) contain many Meta-Analysis objects... '''
    def __init__(self, name, MAs=None):
        self.name = name
        self.MAs = MAs
        
class MetaAnalysis:
    ''' MetaAnalyses (potentially) contain many Study objects... '''
    def __init__(self, name, data_type, method, studies=None):
        self.name = name
        self.data_type = data_type
        self.method = method
        self.studies = studies
        
    def get_method_str(self):
        return methods_d[self.method]

class Study:
    ''' ... This is where we bottom out; the raw data lives here, too. '''
    def __init__(self, index, data_type, raw_data, author, year):
        self.index = index # the index is just the line number 
        self.data_type = data_type
        self.raw_data = raw_data
        self.author = author
        self.year = year
        
        
####################################
#     Here is where we run the     #
#      analyses.                   #
####################################
def run_binary_analysis(studies, method_key):
    # map the method_key (see methods_d for semantics) to the 
    # corresponding function string (including pertinent measure,
    # i.e., metric) in metafor
    metafor_methods_d = {1: "rma.peto(",
                         2: "rma.mh(measure='OR',",
                         3: "rma.uni(measure='OR',",
                         4: "rma.mh(measure='RR',",
                         5: "rma.uni(measure='RR',",
                         6: "rma.mh(measure='RD',",
                         7: "rma.uni(measure='RD',"
    }
    
    '''
    Following metafor's conventions. i.e.,:
    
                outcome 1  outcome 2  total
       group 1    'ai'       'bi'     'n1i'
       group 2    'ci'       'di'     'n2i'
    '''
    ai_str = "c(" + ",".join(\
        ["%s" % study.raw_data[0] for study in studies]) + ")"
    # raw_data[1] is N; we want the `other' events, i.e., N-n1
    bi_str = "c(" + ",".join(\
        ["%s" % (study.raw_data[1] - study.raw_data[0]) for study in studies]) + ")"
    ci_str = "c(" + ",".join(\
        ["%s" % study.raw_data[2] for study in studies]) + ")"
    di_str = "c(" + ",".join(\
        ["%s" % (study.raw_data[3] - study.raw_data[2]) for study in studies]) + ")"
    
    # now stitch the string together.
    r_str = "%s<-%sai=%s,bi=%s,ci=%s,di=%s)" % \
            ("bin_result", metafor_methods_d[method_key], \
             ai_str, bi_str, ci_str, di_str)
    
    #pdb.set_trace()
    return _rls_to_pyd(ro.r(r_str))

def zero_not_in_range((lower, upper)):
    return lower > 0 or upper < 0
    
def run_subsets_for_ma(ma):
    '''
    This is where the magic happens.
    '''
    method_key = ma.method
    N = len(ma.studies)
    # first run a meta-analysis over the original studies
    overall_ma_res = run_binary_analysis(ma.studies, method_key)
    overall_results = parse_results(overall_ma_res, N)
    
    original_was_sig = zero_not_in_range(overall_results["ci"])
    original_q_val_less_than_ten = overall_results["Q"] <= .1
    original_I2_below_fifty = overall_results["I2"] <= 50
    
    absolute_cut_offs = [2.0, 1.5, 1.2]
    exclude_up_to = int(round(.2*N))
    count_ds = {} # map the number excluded to counts
    for n in range(exclude_up_to):
        print "on ma %s, excluding all subsets up to size %s (out of %s)." % (ma.name, n+1, exclude_up_to)
        n_result_d = {"est_stat_sig":0, "Q_stat_sig":0, "I2":0}
        # add entries for the different cut offs
        for cut_off in absolute_cut_offs:
            n_result_d[cut_off] = 0
        subsets = sample_subsets(ma.studies, n+1)
        for subset in subsets:
            cur_results = parse_results(run_binary_analysis(subset, method_key), N-(n+1))
            # now tally up quantities/changes of interest.
            # first, we'll check the point estimate, against our 
            # different cutoffs
            est_diff = abs(overall_results["est"] - cur_results["est"])
            for cut_off in absolute_cut_offs:
                # log the cut off if we're not using risk-difference
                delta = math.log(cut_off) if not "RD" in methods_d[method_key] else cut_off
                if est_diff >= delta:
                    n_result_d[cut_off] += 1
                    
            # next, we'll look at the statistical sig., ie., 
            # whether 0 
            cur_is_sig = zero_not_in_range(cur_results["ci"])
            if cur_is_sig != original_was_sig:
                n_result_d["est_stat_sig"] += 1
                
            # and now we'll look at the Q stat
            cur_q_val_less_than_ten = cur_results["Q"] <= .1
            if cur_q_val_less_than_ten != original_q_val_less_than_ten:
                n_result_d["Q_stat_sig"] += 1
            
            # Finally, have a looksee at I^2
            cur_I2_below_fifty = cur_results["I2"] <= 50
            if cur_I2_below_fifty != original_I2_below_fifty:
                n_result_d["I2"] += 1
                
        count_ds[n+1] = n_result_d
    
    return count_ds    
        
def parse_results(results_d, N):
    d = {}
    d["est"] = results_d["b"]["intrcpt"]
    d["ci"] = (results_d["ci.lb"], results_d["ci.ub"])
    d["Qp"] = results_d["QEp"]
    d["Q"] = results_d["QE"]
    d["I2"] =  max(0,100*((d["Q"]-(N-1))/d["Q"]))
    return d
    
    
def filter(reviews, selection_f):
    ''' 
    Returns a set of meta-analyses from the set of reviews.
    Reviews map to many MAs. Some of these may be correlated
    (i.e., contain some of the same studies) so here we select
    which MAs to keep per each review, as specified by the 
    parametric selection_f.
    '''
    selected_MAs = []
    for review in reviews:
        selected_MAs.append(select_MA(review.MAs))
    
def default_selection(MAs, binary=True):
    '''
    Select a single MA from each review, based on the following algorithm:
    1. If one MA has more studies than all others, select it.
    2. If multiple MAs have the same (max) numbers of studies,
        select from these the one with the most patients
    3. If binary, take the number with largest # of events
    4. Finally, choose at random amongst top contenders
    '''
    candidates = [ma for ma in MAs if ma.studies is not None]
    n_studies_in_candidates = [len(ma.studies) for ma in candidates]
    max_n = max(n_studies_in_candidates)
    candidates = [ma for ma in candidates if len(ma.studies) == max_n]

    if len(candidates)==1:
        return candidates[0]
    # otherwise, two or more MAs have the same (max) num of studies
    # so we look at patient counts
    count_patients = lambda ma: sum([sum(study.raw_data) for study in ma.studies])
    patient_counts_in_candidates = [count_patients(ma) for ma in candidates]
    max_patients = max(patient_counts_in_candidates)
    candidates = [ma for ma in candidates if count_patients(ma) == max_patients]
    # again see if we've whittled it down to 1 study
    if len(candidates)==1:
        return candidates[0]
    # still no dice. if binary, goto step 3.
    if binary:
        count_events = lambda ma: sum([study.raw_data[0] + study.raw_data[2] for study in ma.studies])
        event_counts = [count_events(ma) for ma in candidates]
        max_events = max(event_counts)
        candidates = [ma for ma in candidates if count_events(ma) == max_events]
        if len(candidates)==1:
            return candidates[0]
    # finally, just pick at random
    return random.choice(candidates)
    
    
####################################
#   Enumerating subsets & all      #
#    that jazz.                    #
####################################

# A note here; these methods return the set of
# sets that *exclude* subsets of size k.
def k_subsets_i(n, k):
    '''
    Yield each subset of size k from the set of intergers 0 .. n - 1
    n -- an integer > 0
    k -- an integer > 0
    '''
    # Validate args
    if n < 0:
        raise ValueError('n must be > 0, got n=%d' % n)
    if k < 0:
        raise ValueError('k must be > 0, got k=%d' % k)
    # check base cases
    if k == 0 or n < k:
        yield set()
    elif n == k:
        yield set(range(n))

    else:
        # Use recursive formula based on binomial coeffecients:
        # choose(n, k) = choose(n - 1, k - 1) + choose(n - 1, k)
        for s in k_subsets_i(n - 1, k - 1):
            s.add(n - 1)
            yield s
        for s in k_subsets_i(n - 1, k):
            yield s
    
def sample_subsets(studies, k, sample_size=200000):
    subsets = []
    max_size = 1000000
    if choose(len(studies), k) < max_size:
        # this returns all subsets of x -- notice
        # that in this case the sample_size argument
        # is ignored!
        all_k_sets = [indices for indices in k_subsets_i(len(studies), k)]
        for index_set in all_k_sets:
            subsets.append(set([studies[i] for i in range(len(studies)) if not i in index_set]))
        return subsets
    else:
        # sample from the subset space
        # this is the dumb way.
        already_drawn = []
        n_subsets = 0 # avoid taking length each iteration
        indices = range(len(studies))
        while n_subsets < sample_size:
            candidate_indices = random.sample(indices, k)
            candidate_indices.sort()
            if not candidate_indices in already_drawn:
                already_drawn.append(candidate_indices)
                subsets.append(set([studies[i] for i in range(len(studies)) if not i in candidate_indices]))
                n_subsets+=1
                if n_subsets % 10 == 0:
                    print "%s subsets drawn so far." % n_subsets
        return subsets




####################################
#     Data parsing (wrangling?)    #
####################################
def parse_bin_data():
    '''
    Converts the flat-db of binary reviews/meta-analyses 
    in the tsv file @ binary_data_path into a list
    of Python objects. 
    '''
    bin_db = [x.split("\t") for x in open(binary_data_path).readlines()]
    bin_headers = bin_db[0]
    bin_db = bin_db[1:]
    method_index = bin_headers.index("method")
    method_str_index = bin_headers.index("methods_string\n")
    review_index = bin_headers.index("reviewCode")
    study_index = bin_headers.index("study")
    MA_index = bin_headers.index("sorter")
    
    cur_review = bin_db[0][review_index]
    cur_MA = bin_db[0][MA_index]
    cur_method = int(bin_db[0][method_index])
    # we'll populate these lists as we go.
    reviews = []
    cur_MAs = [] # reviews map to many MAs
    cur_studies = [] # ... which in turn map to studies
    ma_ids_list = []
    # walk over the whole file
    for i, line in enumerate(bin_db):
        if line[MA_index] != cur_MA:
            # instantiate an MA object here, add to the list 
            # of them
            cur_MAs.append(MetaAnalysis(cur_MA, BINARY, cur_method, cur_studies))
            cur_MA = line[MA_index]
            cur_method = int(line[method_index])
            cur_studies = [bin_line_to_study(line, i)]
            
            ## are we on a new review??
            if line[review_index] != cur_review:
                # then we're through with the previous 
                # review; add it to the collection
                reviews.append(Review(cur_review, cur_MAs))
                cur_MAs = [] # this will be appended to...
                cur_review = line[review_index]
        else:
            # otherwise; same MA, same review, so just
            # append the current study to the list
            cur_studies.append(bin_line_to_study(line, i))
    
    reviews.append(Review(cur_review, cur_MAs))
    cur_MAs.append(MetaAnalysis(cur_MA, BINARY, cur_studies))
    return reviews
    
def bin_line_to_study(line, index):
    ### warning -- hard-coding the indices here
    study_au_yr_index = 6
    n1_index, N1_index, n2_index, N2_index = range(8,12)
    raw_data_indices = [n1_index, N1_index, n2_index, N2_index]
    au_yr = line[study_au_yr_index].split(" ")
    study_yr = au_yr[-1]
    study_au = " ".join(au_yr[:-1])
    study_raw_data = [float(line[x]) for x in raw_data_indices]
    return Study(index, BINARY, study_raw_data, study_au, study_yr)
  
    

def _rls_to_pyd(r_ls):
    '''
    This is a hacky way of translating R representations
    from rpy2 to Python native types
    '''
    # base case is that the type is a native python type, rather
    # than an Rvector
    d = {}

    for name, val in zip(r_ls.getnames(), r_ls):
        try:
        
            # first check the key
            if str(name) != "NULL":
                #print name
                #print type(name)
                if "rpy2.robjects" in str(type(name)):
                    name = str(name[0])
                if not "rpy2.robjects" in str(type(val)):
                    # base case; not an rtype object
                    d[name] = val
                elif str(val)=="NULL":
                    d[name] = None
                elif str(val.getnames())=="NULL":
                    d[name] = val[0]
                else:
                    d[name] = _rls_to_pyd(val)
                if not isinstance(name, str):
                    raise Exception, "arg"
            else:
                # name is null
                return val

        except Exception,  inst:
            print "error parsing R tuple.. "
            print inst
            print "ignoring."

    return d
    
def run_binary_analyses(start=0, end=None):
    bin_reviews = parse_bin_data()
    if end is None:
    	end = len(bin_reviews)-1
    selected_MAs = [default_selection(bin_review.MAs) for bin_review in bin_reviews]
    selected_MAs = [ma for ma in selected_MAs if len(ma.studies)>=5]
    output_fields = [1.2, 1.5, 2.0, 'I2', 'est_stat_sig', 'Q_stat_sig']
    fout = open("binary_results_%s_%s" % (start, end), 'w')
    fout.write("\t".join(["MA_ID", "k", "nstudies", "method_str", \
                    "1.2", "1.5", "2.0", 'I2', 'est_stat_sig', 'Q_stat_sig']) + "\n")
    for ma in selected_MAs[start:end]:
        res_d = run_subsets_for_ma(ma)
        ns = res_d.keys()
        ns.sort()
        for k in ns:
            cur_line = [ma.name, k, len(ma.studies), methods_d[ma.method]]
            cur_counts = res_d[k]
            cur_line.extend([cur_counts[1.2], cur_counts[1.5], cur_counts[2.0], \
                             cur_counts["I2"], cur_counts["est_stat_sig"], cur_counts["Q_stat_sig"]])
            cur_line = "\t".join([str(x) for x in cur_line])+"\n"                    
            fout.write(cur_line)
    fout.close()
    
def choose(n, k):
    if 0 <= k <= n:
        p = 1
        for t in xrange(min(k, n - k)):
            p = (p * (n - t)) // (t + 1)
        return p
    else:
        return 0

if __name__ == '__main__':
	start = int(sys.argv[1])
	end = int(sys.argv[2])
	run_binary_analyses(start, end)
	

            