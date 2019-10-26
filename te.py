from ortools.linear_solver import pywraplp
from random import randint
import pandas

K = [
    {
        "fields" : [0,1,2,3],
        "addr"   : [],
        "search" : [],
    },
    {
        "fields" : [4,5,6,7],
        "addr"   : [],
        "search" : [],
    },
    {
        "fields" : [8,9,10,11],
        "addr"   : [],
        "search" : [],
    },
    {
        "fields" : [12,13,14,15],
        "addr"   : [],
        "search" : [],
    },
]

F = [128 for _ in range(16)]

def SolVal(x):
  if type(x) is not list:
    return 0 if x is None \
      else x if isinstance(x,(int,float)) \
           else x.SolutionValue() if x.Integer() is False \
                else int(x.SolutionValue())
  elif type(x) is list:
    return [SolVal(e) for e in x]

def solve_model(KeyInfo, field_size):
    s    = pywraplp.Solver("KeyMaker", pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)
    #------------------------------------------------------------
    #        *key field usage
    #             key_field_mbr[key][field] = [0,1] (Input)
    #------------------------------------------------------------
    numKey   = len(KeyInfo)
    numField = len(field_size)
    numKeyMakerUnit = 4
    numPHVFlit      = 12
    numProfile      = 8
    keyMakerSize    = 256
    maxKeySize      = 512
    #------------------------------------------------------------
    #        *key field size check
    #             sum(key_field_mbr[key][field] * field_size[field] for all field) < 512 (input check)
    #------------------------------------------------------------
    for key in range(numKey):
        total_size = 0
        for field in KeyInfo[key]['fields']:
            total_size += field_size[field]
        if(total_size > 512):
            print("ERROR: key " + str(key) + " size > "+ str(maxKeySize) +" bits not supported\n")
            return
    #------------------------------------------------------------
    #        *key maker hardware unit field membership
    #             km_hw[flit][unit][field] = [0,1]  (V)
    #------------------------------------------------------------
    km_hw = [[[s.IntVar(0,1,'') for f in range(numField)] for unit in range(numKeyMakerUnit)] for flit in range(numPHVFlit)]
    #------------------------------------------------------------
    #        *key maker hardware unit capacity
    #             sum(km_hw[flit][unit][field] * field_size[field]) < 256  (C)
    #------------------------------------------------------------
    for flit in range(numPHVFlit):
        for unit in range(numKeyMakerUnit):
            s.Add(sum(km_hw[flit][unit][field] for field in range(numField)) <= keyMakerSize - 1)
    #------------------------------------------------------------
    #        *assign a field to a phv flit
    #             sum(km_hw[flit][unit][field] <= 1 for all flit and unit (C)
    #------------------------------------------------------------
    phv = [[s.IntVar(0,1,'') for f in range(numField)] for flit in range(numPHVFlit)]
    for field in range(numField):
        s.Add(sum(phv[flit][field] for flit in range(numPHVFlit)) == 1)
    #------------------------------------------------------------
    #        *key maker update
    #             phv[flit][field] = [0,1] (V)
    #             km_hw[flit][unit][field] <= km_hw[flit-1][unit][field]) + phv[flit][field] (C)
    #------------------------------------------------------------
    for flit in range(numPHVFlit):
        for unit in range(numKeyMakerUnit):
            for field in range(numField):
                if flit == 0 :
                    s.Add(km_hw[flit][unit][field] <= phv[flit][field])
                else :
                    s.Add(km_hw[flit][unit][field] <= km_hw[flit-1][unit][field] + phv[flit][field])
    #------------------------------------------------------------
    #        *key maker hw unit selection
    #             km_hw_sel[key][flit][unit] = [0,1] (V)
    #             sum(km_hw_sel[key][flit][unit] for all unit) >= 1 for all key (C)
    #             sum(km_hw_sel[key][flit][unit] for all unit) <= 2 for all key (C)
    #------------------------------------------------------------
    km_hw_sel = [[[s.IntVar(0,1,'') for unit in range(numKeyMakerUnit)] for flit in range(numPHVFlit)] for k in range(numKey)]
    for key in range(numKey):
        s.Add(sum(km_hw_sel[key][flit][unit] for flit in range(numPHVFlit) for unit in range(numKeyMakerUnit)) >= 1)
        s.Add(sum(km_hw_sel[key][flit][unit] for flit in range(numPHVFlit) for unit in range(numKeyMakerUnit)) <= (maxKeySize/keyMakerSize))
        for unit in range(numKeyMakerUnit):
            s.Add(sum(km_hw_sel[key][flit][unit] for flit in range(numPHVFlit)) <= 1)
    for flit in range(numPHVFlit):
        for unit in range(numKeyMakerUnit):
            s.Add(sum(km_hw_sel[key][flit][unit] for key in range(numKey)) <= 1)

    #------------------------------------------------------------
    #        *key profile usage (assume single byte mode for now)
    #             km_hw_prof[flit][unit] = [0,1] (V)
    #             sum(km_hw_prof[flit][unit] for all flit/unit) <= 8 (C)
    #------------------------------------------------------------
    km_hw_profile = [[s.IntVar(0,1,'') for unit in range(numKeyMakerUnit)] for flit in range(numPHVFlit)]
    s.Add(sum(km_hw_profile[flit][unit] for flit in range(numPHVFlit) for unit in range(numKeyMakerUnit)) <= numProfile - 1)

    #------------------------------------------------------------
    #        *key maker can only select hw unit that has a profile
    #             km_hw_sel[key][flit][unit] <= km_hw_prof[flit][unit] for all key,flit,unit (C)
    #------------------------------------------------------------
    for key in range(numKey):
        for flit in range(numPHVFlit):
            for unit in range(numKeyMakerUnit):
                s.Add(km_hw_sel[key][flit][unit] <= km_hw_profile[flit][unit])
    #------------------------------------------------------------
    #        *key maker has the fields required
    #             km_hw[flit][unit][field] >= km_hw_sel[key][flit][unit] for all key, flit, unit, key_field_mbr[key][field]
    #------------------------------------------------------------
    for key in range(numKey):
        for field in KeyInfo[key]['fields']:
            for flit in range(numPHVFlit):
                for unit in range(numKeyMakerUnit):
                    s.Add(km_hw[key][unit][field] >= km_hw_sel[key][flit][unit])
    #------------------------------------------------------------
    #        *key field assignment
    #             km_field_byte_loc[key][field] = [0,31]  (V)
    #             km_field_bit_loc[key][field]  = [0,7]  (V)
    #             (km_field_byte_loc[key][field] * 8 + km_field_bit_loc[key][field]) + field_size[field] < 511  (C) for all key,field
    #------------------------------------------------------------
    km_field_byte_loc = [[s.IntVar(0,31,'') for field in range(numField)] for key in range(numKey)]
    km_field_bit_loc  = [[s.IntVar(0,7,'')  for field in range(numField)] for key in range(numKey)]
    for key in range(numKey):
        for field in KeyInfo[key]['fields']:
            s.Add((km_field_byte_loc[key][field] * 8) + km_field_bit_loc[key][field] + field_size[field] <= maxKeySize - 1)
    #------------------------------------------------------------
    #        *key address field needs to be 16 bit aligned
    #             km_field_addr_sel[key][field] = [0,31]  (V)
    #             km_field_start[key][field] + field_size[field] = key_field_addr_sel[key_field] * 16 (C)
    #------------------------------------------------------------
    km_field_aadr_sel = [[s.IntVar(0,31,'') for field in range(numField)] for key in range(numKey)]
    for key in range(numKey):
        for field in KeyInfo[key]['addr']:
            s.Add((km_field_byte_loc[key][field] * 8) + km_field_bit_loc[key][field] + field_size[field] == km_field_aadr_sel[key][field] * 16)
    #------------------------------------------------------------
    #        *key search field needs to be grouped together
    #             search_size = sum(km_field_search[key][field])
    #             abs(km_field_search[key][field_x],km_field_search[key][field_y]) < search_size for all fields)
    #                         or
    #             km_field_search[key][field_x] - km_field_search[key][field_y]) < search_size (C)
    #             km_field_search[key][field_y] - km_field_search[key][field_x]) < search_size (C)
    #------------------------------------------------------------
    for key in range(numKey):
        search_size = sum(field_size[KeyInfo[key]['search'][field]] for field in KeyInfo[key]['search'])
        for fieldx in KeyInfo[key]['search']:
            for fieldy in KeyInfo[key]['search']:
                if (fieldx != fieldy):
                    s.Add(((km_field_byte_loc[key][fieldx] * 8) + km_field_bit_loc[key][fieldx] + field_size[fieldx] -1) - ((km_field_byte_loc[key][fieldy] * 8) + km_field_bit_loc[key][fieldy]) <= search_size - 1)
                    s.Add(((km_field_byte_loc[key][fieldy] * 8) + km_field_bit_loc[key][fieldy] + field_size[fieldy] -1) - ((km_field_byte_loc[key][fieldx] * 8) + km_field_bit_loc[key][fieldx]) <= search_size - 1)
    #------------------------------------------------------------
    #        *key maker bit select usage
    #             km_field_bitsel0[key][field][bit] = [0,1] (V)
    #             km_field_bitsel1[key][field][bit] = [0,1] (V)
    #             sum(km_field_bitsel0[key][field][bit]) < 8 (C)
    #             sum(km_field_bitsel1[key][field][bit]) < 8 (C)
    #------------------------------------------------------------
    km_field_bitsel0  = [[[s.IntVar(0,1,'')  for bit in range(field_size[field])] for field in range(numField)] for key in range(numKey)]
    km_field_bitsel1  = [[[s.IntVar(0,1,'')  for bit in range(field_size[field])] for field in range(numField)] for key in range(numKey)]
    s.Add(sum(km_field_bitsel0[key][field][bit] for key in range(numKey) for field in range(numField) for bit in range(field_size[field]) ) <= 7)
    s.Add(sum(km_field_bitsel1[key][field][bit] for key in range(numKey) for field in range(numField) for bit in range(field_size[field]) ) <= 7)
    #------------------------------------------------------------
    #        *key field bits alignments (keep them consecutve for now)
    #             km_field_align[key][field][bit] = [0,7] (V)
    #             km_field_byte_wrap[key][field][bit] = [0,field_size/8] (V) (for modolo operation)
    #             for (i=1; i < field_size[field]; i++)
    #                km_field_align[key][field][i] = (km_field_bit_loc[key][field] + i) - (8 * km_field_byte_wrap[key][field][bit])
    #------------------------------------------------------------
    km_field_align = {}
    km_field_wrap = {}
    for key in range(numKey):
        km_field_align[key] = {}
        km_field_wrap[key] = {}
        for field in KeyInfo[key]['fields']:
            km_field_align[key][field] = {}
            km_field_wrap[key][field] = {}
            for bit in range(field_size[field]):
                km_field_align[key][field][bit] = s.IntVar(0,7,'')
                km_field_wrap[key][field][bit]  = s.IntVar(0,field_size[field]/8,'')
                s.Add(km_field_align[key][field][bit] == (km_field_bit_loc[key][field] + bit) - (km_field_wrap[key][field][bit] * 8))
    #------------------------------------------------------------
    #        *key maker byte align with phv
    #             phv_field_align[field][bit]     = [0,7] (V)
    #             km_field_align[key][field][bit] == phv_field_align[field][bit] if km_field_bitsel0[key][field][bit] == 0 &&  km_field_bitsel1[key][field][bit] == 0
    #                        or
    #             km_field_align[key][field][bit] - phv_field_align[field][bit]      <= km_field_bitsel0[key][field][bit] * X + km_field_bitsel1[key][field][bit] * X (C)
    #             phv_field_align[field][bit]     - km_field_align[key][field][bit]  <= km_field_bitsel0[key][field][bit] * X + km_field_bitsel1[key][field][bit] * X (C)
    #                  note: * X is an arbitrary big number (e.g. 8)
    #                        * if any bit select is used, want this equation to always be true (i.e. constraint not applied) 
    #------------------------------------------------------------
    phv_field_align = {}
    for field in range(numField):
        phv_field_align[field] = {}
        for bit in range(field_size[field]):
            phv_field_align[field][bit] = s.IntVar(0,7,'')
    for key in range(numKey):
        for field in KeyInfo[key]['fields']:
            for bit in range(field_size[field]):
                #print("key=%d, field=%d, bit=%d\n" % (key,field,bit) )
                s.Add(km_field_align[key][field][bit] - phv_field_align[field][bit] <= (km_field_bitsel0[key][field][bit] * 8) + (km_field_bitsel1[key][field][bit] * 8))
    rc = s.Solve()

    print("rc=", rc)
    #print("km_hw=", rc)
    #print(pandas.DataFrame(SolVal(km_hw)))
    for key in range(numKey):
        print("km_hw_sel: key=", key)
        print(pandas.DataFrame(SolVal(km_hw_sel[key])).T.replace({0 : '.'}))
        print("")
    print("phv field location")
    print(pandas.DataFrame(SolVal(phv)).T.replace({0 : '.'}))
    print("")
    #print("phv field alignment")
    #print(pandas.DataFrame(SolVal(phv_field_align)).T.replace({0 : '.'}))
    #print("")
    print("key maker hardware profile")
    print(pandas.DataFrame(SolVal(km_hw_profile)).T.replace({0 : '.'}))
    print("")
    print("key maker byte location")
    for key in range(numKey):
        print(pandas.DataFrame(SolVal(km_field_byte_loc)))
        print("")
    print("key maker bit location")
    for key in range(numKey):
        print(pandas.DataFrame(SolVal(km_field_bit_loc)))
        print("")


    print("km_hw=")
    print(pandas.DataFrame(SolVal(km_hw)).replace({0 : '.'}))
    print("")
    
#    rnb = SolVal(nb)
#    return rc,rnb,rolls(rnb,SolVal(x),SolVal(w),D),SolVal(w)


solve_model(K,F)


    
