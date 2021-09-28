
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 14 19:10:45 2019

@author: jmw
"""
############################################
import numpy as np
import matplotlib.pyplot as plt
import csv
import scipy.stats as st
from matplotlib.lines import Line2D
import seaborn as sns
import pandas as pd
from statannot import add_stat_annotation
import matplotlib
matplotlib.rc('font', family='sans-serif') 
matplotlib.rc('font', serif='Arial')
##############################################
name_score = 'AAV1_PG_300 trials.csv' ## user scored epochs
days = 2 ## number of days per mouse
bout_length = 5 ## number of seconds per bout
globalFont = 15 ## fontsize for plotting
global_name = name_score[0:len(name_score)-4]  ## name attached to 
min_before1st = 3 #min before the first trial
beforeAndAfter = 3 #min before and after to visualize
num_trials = 300 
inhibitory = False ## if True the lights on is yellow, otherwise blue
ITERATIONS = 10000 ## number of iterations for bootstrap
bout_per_min = int(60/bout_length)
stim_freq = 20 ## how ofter stim repeats in min
stim_length = 90 ##stim length in sec
trials_per_animal = 50 ## not really used now that it calculates automatically

def extract_column(a, column):
    #take a single column from a more complex matrix
    ## returns that column
    new_matrix = []
    for i in range(len(a)):
        new_matrix.append(a[i][column])
    return new_matrix 

def matrix_float(a):
    #converts a matrix of strings to floats
    ##often necessary when importing the data csv files into python
    new_matrix = []
    for i in range(len(a)):
        try:
            new_matrix.append((float(a[i])))
        except ValueError:
            print(i)
            print("matrixFloat")
    return new_matrix

def extract_row(a, column):
    #take a single column from a more complex matrix
    ## returns that column
    new_matrix = []
    for i in range(len(a[0])):
        new_matrix.append(a[column][i])
    return new_matrix 

def createDataFrame(frame, data, name, index):
    ## creates pd dataframe with the given parameters
    for i in range(len(data)):
        frame[name][i+index] = data[i]
    return frame

def plotDots_ttest(wt,wtName, mut,mutName, ylabel,name):
    ## plots a boxplot and swarmplot for 2 groups    
    labels = []
    for i in range(len(wt)):
        labels.append(wtName)
    for i in range(len(mut)):
        labels.append(mutName)
    height = len(wt) + len(mut)
    framee = pd.DataFrame(np.random.randint(low=0.0, high= 100, size=(height, 2)), columns=['cFos', 'condition'],dtype = np.float64)
    framee = createDataFrame(framee, wt, 'cFos', 0)
    framee = createDataFrame(framee, mut, 'cFos', len(wt))
    framee = createDataFrame(framee, labels, 'condition',0)
    sns.set_style("ticks")
    sns.set(font_scale = 1.4)
    sns.set_style("ticks")
    g2 = sns.boxplot(x="condition", y="cFos", data=framee, palette = ["khaki",'mediumorchid'], linewidth =1)
    sns.swarmplot(x = "condition", y = "cFos", data = framee,palette = ["goldenrod",'blueviolet'], alpha = 1, s = 8)
    sns.despine()
    order = [wtName, mutName]
    test_results = add_stat_annotation(g2, data=framee, x='condition', y='cFos', order=order,
                                   box_pairs=[(wtName, mutName)],
                                   test='t-test_ind', text_format='star',
                                   loc='outside', verbose=2)
    plt.xlabel('')
    plt.ylabel(ylabel, fontsize = globalFont)
    total_wake_ttest = st.ttest_ind(wt,mut, equal_var = False)
    total_wake_ttest = total_wake_ttest[1]
    # plt.title(str(total_wake_ttest), fontsize = globalFont)
    plt.savefig(name +global_name+".pdf")
    plt.show()
    return

def scrub_formating(a):
    ## gets to right formatting from the extracted data
    stop = len(a[0])-1
    for i in range(len(a)):
        for j in range(len(a[0])):
            
            a[i][j] = a[i][j][2:len(a[i][j])-1]
    for i in range(len(new_data)):
        new_data[i][stop] = new_data[i][stop][0:len(new_data[i][stop])-1]
    a[0][0] = a[0][0][1:len(a[0][0])] 
    a[0][stop] = a[0][stop][0:len(a[0][stop])-1]     
    return a

def just_scores(a):
    ## gets just the scores from the raw data
    for i in range(len(a)):
            a[i] = a[i][2:len(a[i])]
    return a

def calc_perc(a, name):
    ## derives percent in a state in a given minute
    counter = 0
    for i in range(len(a)):
        if a[i] == name:
            counter +=1
    perc = counter/len(a)*100
    return perc

def calc_percent_tot(a):
    ## cacluates the percent for every bout timepoint
    wake_perc = []
    NREM_perc= []
    REM_perc = []
    new_matrix = []
    if len(a) ==1:
        for i in range(len(NREM_avg)):
            new_matrix.append(0)
        return new_matrix, new_matrix, new_matrix
    for i in range(len(a[0])):
        column = extract_column(a, i)
        curr_wake = calc_perc(column, 'Wake')
        curr_NREM = calc_perc(column, 'NREM')
        curr_REM = calc_perc(column, 'REM')
        wake_perc.append(curr_wake)
        NREM_perc.append(curr_NREM)
        REM_perc.append(curr_REM)
    return wake_perc, NREM_perc, REM_perc

def create_sec():
    ## create seconds for the plot
    tot_bout = int(2*bout_per_min*beforeAndAfter + stim_length/bout_length)
    sec = []
    for i in range(tot_bout):
        sec.append(bout_length*i)
    for i in range(len(sec)):
        sec[i] = sec[i] - bout_per_min*beforeAndAfter*bout_length
    return sec


def create_sec_CI():
    ## create seconds for the plot
    tot_bout = int(2*bout_per_min*beforeAndAfter + stim_length/bout_length)
    sec = []
    for i in range(tot_bout):
        sec.append(bout_length*i)
    for i in range(len(sec)):
        sec[i] = sec[i] - bout_per_min*beforeAndAfter*bout_length + 2*bout_length
    return sec

def plotPercSEM(wake, wakeSEM, NREM,NREMSEM,  REM,REMSEM, sec, name):
    ## plots percent wake, NREM, REM with the standard error
    ax1 = plt.axes(frameon=False)
    plt.errorbar(sec, wake,wakeSEM, c = 'blue', lw = 4, ecolor= 'dodgerblue')
    plt.errorbar(sec, NREM,NREMSEM, c = 'purple', lw = 4, ecolor= 'thistle')
    plt.errorbar(sec, REM,REMSEM, c = 'gray', lw = 4, ecolor= 'lightgrey')
    plt.xlabel('Time (min)',fontsize = globalFont)
    plt.ylabel('%',fontsize = globalFont)    
    plt.ylim(0,100)
    xmin, xmax = ax1.get_xaxis().get_view_interval()
    ymin, ymax = ax1.get_yaxis().get_view_interval()
    ax1.tick_params(labelsize=globalFont*.8)  
    plt.xticks(ticks = [-180,-120,-60, 0, 60, 120, 180, 240],labels = ['-3','-2','-1','0','1','2','3', '4'])
    ax1.add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=2))
    ax1.add_artist(Line2D((xmin, xmin), (ymin, ymax), color='black', linewidth=2))  
    if inhibitory == False:
        ax1.axvspan(0,stim_length, color = 'lightskyblue')
    if inhibitory == True:
         ax1.axvspan(0,stim_length, color = 'yellow')       #plt.title(name)
    plt.ylim(0,102)
    #plt.title(name)
    plt.savefig(name +'opto_stim.pdf')
    plt.show()
    return

def extract_single_stage_before(a, name):
    ## extract just a certain sleep stage before
    ## name is the state like 'REM', 'NREM', or 'Wake'
    matrix = []
    index = int(bout_per_min*beforeAndAfter)
    for i in range(len(a)):
        if a[i][index] == name:
            matrix.append(a[i])
    new_matrix = []
    if len(matrix) == 0:
        new_matrix.append(0)
        matrix = new_matrix
    return matrix

def perc_per_animal_variable(a, TRIALS):
    ## finds the percent per animal if there is a variable number of trials per animal
    wake_tot = []
    NREM_tot = []
    REM_tot = []
    num_animals = int(len(TRIALS))
    start = 0
    stop = 0
    for i in range(num_animals):
        start = stop
        stop = stop + TRIALS[i]
        WAKEPER, NREMPER, REMPER = calc_percent_tot(a[start:stop])
        wake_tot.append(WAKEPER)
        NREM_tot.append(NREMPER)
        REM_tot.append(REMPER)
    return wake_tot, NREM_tot, REM_tot

def trials_by_animal(a, TRIALS):
    ## calculates the trials per animal and returns a trial list
    b = []
    start = 0
    stop = 0
    num_animals = int(len(TRIALS))
    for i in range(num_animals):
        start = stop
        stop = stop + TRIALS[i]
        curr_b = a[start:stop]
        b.append(curr_b)
    return b
    
def cleanup_neg_trials(a):
    ## deletes zeros from a matrix
    count = 0
    for i in range(len(a)):
        all_zeros = np.sum(a[i-count])
        if all_zeros == 0:
            curr_delete = i-count
            a.pop(curr_delete)
            count+=1
    return a

def perc_per_animal_single_stage_variable(a, TRIAL, name):
    ## finds the percent per stage when there is a variable # of trials per mouse
    wake_tot = []
    NREM_tot = []
    REM_tot = []
    num_animals = int(len(TRIAL))
    total_trials = 0
    start = 0
    stop = 0
    for i in range(num_animals):
        start = stop
        stop = stop + TRIAL[i]
        curr_data = a[start:stop]
        curr_before = extract_single_stage_before(curr_data, name)
        total_trials = total_trials + len(curr_before)
        WAKEPER, NREMPER, REMPER = calc_percent_tot(curr_before)
        wake_tot.append(WAKEPER)
        NREM_tot.append(NREMPER)
        REM_tot.append(REMPER)
    wake_tot = cleanup_neg_trials(wake_tot)
    NREM_tot = cleanup_neg_trials(NREM_tot)
    REM_tot = cleanup_neg_trials(REM_tot)
    return wake_tot, NREM_tot, REM_tot, total_trials

def avg_animals(a):
    ## basically for a given nested matrix of multiple trials
    ##returns 2 single matrixes: the average and one of standard deviation
    ## there's probably a python function that does this already but made my own 
    avg_matrix = []
    std_matrix = []
    if len(a) == 0:
        print('this metric has zero occurances')
        avg_matrix = 0
        for i in range(len(avg_matrix)):
            avg_matrix[i] = 0
        std_matrix = avg_matrix
        return avg_matrix, std_matrix
    for i in range(len(a[0])):
        curr_total = []
        for j in range(len(a)):
            curr_total.append(a[j][i])
        avg_matrix.append(np.mean(curr_total))
        std_matrix.append(st.sem(curr_total))
    return avg_matrix, std_matrix      

def determineTrialsPerAnimal(trials):
    ## basically 
    diff_trials = []
    counter = []
    for i in range(len(trials)):
        if trials[i][0] == "'":
            trials[i] = trials[i][1:len(trials[i])]
    for i in range(len(trials)):
        if trials[i] not in diff_trials:
            diff_trials.append(trials[i])
    trials = np.asarray(trials)
    for i in range(len(diff_trials)):
        new_num = (trials == diff_trials[i]).sum()
        counter.append(new_num)
    return counter

def undueBootstrapFormatting(a):
    ## undoes the bootstrap formatting in order to do other analysis
    b = []
    for i in range(len(a)):
        for j in range(len(a[i])):
            b.append(a[i][j])
    return b

def all_mouse_single_bootstrap_better(a):
    ## randomly draws out trials with replacement
    num_trials = len(a)
    b = []
    for i in range(num_trials):
        index = int(np.random.rand()*num_trials)
        b.append(a[index])
    return b

def state_probabilities_single_bootstrap_better(a, trialOcc):
    ## for a single bootsrap, gets the random draw
    ## then returns the state probabilities 
    trials4bootstrap = trials_by_animal(a,trialOcc)
    trials4bootstrap = undueBootstrapFormatting(trials4bootstrap)
    curr_bootstrap = all_mouse_single_bootstrap_better(trials4bootstrap)
    wake, NREM, REM = perc_per_animal_variable(curr_bootstrap, trialOcc)
    wake_tot, NREM_tot, REM_tot =  calc_percent_tot(curr_bootstrap)
    return wake_tot, NREM_tot, REM_tot

def combine2matrix(a,b):
    ## combines 2 matrixes into 1
    c = []
    for i in range(len(a)):
        c.append(a[i])
    for i in range(len(b)):
        c.append(b[i])
    return c

def combineIterations(a, b):
    ## combines the iterations from bootstrapping into the correct formatting
    c = []
    for i in range(len(b)):
        if len(a) >0:
            new = combine2matrix(a[i],b[i])
            c.append(new)
        else: 
            return b
    return c

def state_probabilities_all_bootstrap_better(a,trialOcc, iterations = ITERATIONS):
    ## bootstrap for all the trials
    wake_all = []
    NREM_all = []
    REM_all = []
    for i in range(iterations):
        curr_wake, curr_NREM, curr_REM = state_probabilities_single_bootstrap_better(a, trialOcc)
        wake_all.append(curr_wake)
        NREM_all.append(curr_NREM)
        REM_all.append(curr_REM)
    return wake_all, NREM_all, REM_all    

def extractConfidenceIntervals_single(a, index, lower = 0.025, upper = 0.975):
    ## determine the 95% confidence interval for a single trial  
    ## to change the interval, adjust the upper and lower 
    curr_a = extract_column(a, index)
    a_sorted = curr_a[:]
    a_sorted.sort()
    lower_val = a_sorted[int(len(a_sorted)*lower)]
    upper_val = a_sorted[int(len(a_sorted)*upper)]
    return lower_val, upper_val

def extractConfidenceIntervals_all(a, lower = 0.025, upper = 0.975):
    ## determine the 95% confidence interval. 
    ## to change the interval, adjust the upper and lower 
    lower_all =[]
    upper_all = []
    for i in range(len(a[0])):
        low_curr, up_curr = extractConfidenceIntervals_single(a,i)
        lower_all.append(float(low_curr))
        upper_all.append(float(up_curr))
    lower_all = np.array(lower_all)
    upper_all = np.array(upper_all)
    CI = [lower_all, upper_all]
    return CI
        
def plotPercCI(wake, wakeCI, NREM ,NREMCI,  REM,REMCI, sec, name):
    ## plot the wake, NREM, REM % along with the confidence interval from the bootstrap
    sec = create_sec_CI()
    ax1 = plt.axes(frameon=False)
    plt.errorbar(sec, wake,yerr= wakeCI, c = 'blue', lw = 2, ecolor= 'dodgerblue',elinewidth =4)
    plt.errorbar(sec, NREM,NREMCI, c = 'purple', lw = 2, ecolor= 'thistle',elinewidth = 4)
    plt.errorbar(sec, REM,REMCI, c = 'gray', lw = 2, ecolor= 'lightgrey', elinewidth = 4)
    plt.xlabel('Time (min)',fontsize = globalFont)
    plt.ylabel('%',fontsize = globalFont)    
    plt.ylim(0,100)
    xmin, xmax = ax1.get_xaxis().get_view_interval()
    ymin, ymax = ax1.get_yaxis().get_view_interval()
    ax1.tick_params(labelsize=globalFont*.8)  
    plt.xticks(ticks = [-180,-120,-60, 0, 60, 120, 180, 240],labels = ['-3','-2','-1','0','1','2','3', '4'])
    ax1.add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=2))
    ax1.add_artist(Line2D((xmin, xmin), (ymin, ymax), color='black', linewidth=2))  
    if inhibitory == False:
        ax1.axvspan(0,stim_length, color = 'lightskyblue')
    if inhibitory == True:
         ax1.axvspan(0,stim_length, color = 'yellow')       #plt.title(name)
    plt.ylim(0,102)
    plt.savefig(name +'opto_stim.pdf')
    plt.show()
    return

def generatePlotableCI(a, err):
    ## takes the confidence interval previous calculated it 
    ## and puts it in the format the plt.plot() will recognize
    for i in range(len(a)):
        err[0][i] = a[i] - err[0][i]
        err[1][i] =err[1][i] - a[i]
    return err

def wakeMaintenance_single(wake,NREM, REM):
    ## for a single mouse, finds the state % before and after stim. In this function -3 to -2min
    ## and +1.5 to +4.5min
    wakeb4 = []
    NREMb4 = []
    REMb4 = []
    wake_stim = []
    NREM_stim = []
    REM_stim = []
    ###find start of stim
    boutPerMinute = int(60/bout_length)
    stim_start = int(beforeAndAfter*60/bout_length+stim_length/bout_length +2*boutPerMinute)
    stim_bout_num = boutPerMinute
    before_start = 0
    for i in range(stim_bout_num):
        wakeb4.append(wake[before_start +i])
        NREMb4.append(NREM[before_start +i])
        REMb4.append(REM[before_start +i])
    for i in range(stim_bout_num):
        wake_stim.append(wake[stim_start +i])
        NREM_stim.append(NREM[stim_start +i])
        REM_stim.append(REM[stim_start +i]) 
    wakeb4 = np.mean(wakeb4)
    NREMb4 = np.mean(NREMb4)
    REMb4 = np.mean(REMb4)
    wake_stim = np.mean(wake_stim)
    NREM_stim = np.mean(NREM_stim)
    REM_stim = np.mean(REM_stim)
    return wakeb4, wake_stim, NREMb4, NREM_stim, REMb4, REM_stim

def wakeMaintenance_all(wake, NREM, REM):
    ## for all mice,  finds the % of wake, NREM, REM 3 minutes before and after stim
    ## then averages and plots it
    wakeb4 = []
    NREMb4 = []
    REMb4 = []
    wake_after = []
    NREM_after = []
    REM_after = [] 
    for i in range(len(wake)):
        WB, WA, NB,NA, RB, RA = wakeMaintenance_single(wake[i], NREM[i], REM[i])
        wakeb4.append(WB)
        NREMb4.append(NB)
        REMb4.append(RB)      
        wake_after.append(WA)
        NREM_after.append(NA)
        REM_after.append(RA)    
    plotDots_ttest(wakeb4, '3 min before stim', wake_after, '3 min after stim', '% Wake', 'first min compared to last min_wake')
    plotDots_ttest(NREMb4, 'first min', NREM_after, 'last min', '% NREM', 'first min compared to last min_NREM')
    plotDots_ttest(REMb4, 'first min', REM_after, 'last min', '% REM', 'first min compared to last min_REM')
    return wakeb4, wake_after

def mean_probabilites(wake, NREM, REM):
    ## Calculate the probability of wake in the pre-stim period and 
    ## during stim. Then subtract the stim period from the pre-stim
    wakeb4 = []
    NREMb4 = []
    REMb4 = []
    wake_stim = []
    NREM_stim = []
    REM_stim = []
    ###find start of stim
    stim_start = int(beforeAndAfter*60/bout_length)
    stim_bout_num = int(stim_length/bout_length)
    before_start = stim_start - stim_bout_num
    for i in range(stim_bout_num):
        wakeb4.append(wake[before_start +i])
        NREMb4.append(NREM[before_start +i])
        REMb4.append(REM[before_start +i])
    for i in range(stim_bout_num):
        wake_stim.append(wake[stim_start +i])
        NREM_stim.append(NREM[stim_start +i])
        REM_stim.append(REM[stim_start +i]) 
    wakeb4 = np.mean(wakeb4)
    NREMb4 = np.mean(NREMb4)
    REMb4 = np.mean(REMb4)
    wake_stim = np.mean(wake_stim)
    NREM_stim = np.mean(NREM_stim)
    REM_stim = np.mean(REM_stim)
    wake_diff = wake_stim-wakeb4
    NREM_diff = NREM_stim-NREMb4
    REM_diff = REM_stim-REMb4
    return wake_diff, NREM_diff, REM_diff

def mean_prob_bootstrap(wake, NREM, REM):
    wake_all = []
    NREM_all = []
    REM_all = []
    for i in range(len(wake)):
        wake_diff, NREM_diff, REM_diff = mean_probabilites(wake[i], NREM[i], REM[i])
        wake_all.append(wake_diff)
        NREM_all.append(NREM_diff)
        REM_all.append(REM_diff)
    return wake_all, NREM_all, REM_all
        
def increaseORdecrease(a):
    ## finds how many + or - values there are in a matrix
    inc = []
    dec = []
    for i in range(len(a)):
        if a[i] > 0:
            inc.append(a[i])
        elif a[i] < 0:
            dec.append(a[i])
    return inc, dec 

 
def add_up(w,n,r):
    ALL = []
    for i in range(len(w)):
        ALL.append(w[i]+n[i]+r[i])
    return ALL

with open(name_score, newline='\n' ) as inputfile:
   data = list(csv.reader(inputfile)) 
   data[0][0] = data[0][0][1:len(data[0][0])]
NREM_before = '_' +global_name + ' NREM before' 
REM_before = '_' +global_name +'REM before'
wake_before = '_' +global_name +'Wake before'
for i in range(len(data)-1):
    data[i+1][0] = data[i+1][0].split(",")
new_data = [0 for i in range(len(data))]
new_data[0] = data[0]
for i in range(len(data)-1):
    new_data[i+1] = data[i+1][0]
new_data = scrub_formating(new_data)
new_data[0][0]
names = extract_column(new_data, 0)
names[0] = names[0][0:len(names[0])]
trialOccurances = determineTrialsPerAnimal(names)
justData = just_scores(new_data)
wake_total, NREM_total, REM_total = perc_per_animal_variable(justData, trialOccurances)
wake_avg, wake_sem = avg_animals(wake_total)
NREM_avg, NREM_sem = avg_animals(NREM_total)
REM_avg, REM_sem = avg_animals(REM_total)
NREMbefore = extract_single_stage_before(justData, 'NREM')
REMbefore = extract_single_stage_before(justData, 'REM')
wakebefore = extract_single_stage_before(justData, 'Wake')
wakeb4, wakeAfter = wakeMaintenance_all(wake_total, NREM_total, REM_total)
sec = create_sec()


ALL = add_up(wake_avg, NREM_avg, REM_avg)

NREMbefore_wake_tot, NREMbefore_NREM_tot, NREMbefore_REM_tot, NREMbefore_totaltrials = perc_per_animal_single_stage_variable(justData, trialOccurances,'NREM')
NREMbefore_wake_avg, NREMbefore_wake_sem = avg_animals(NREMbefore_wake_tot)
NREMbefore_NREM_avg, NREMbefore_NREM_sem = avg_animals(NREMbefore_NREM_tot)
NREMbefore_REM_avg, NREMbefore_REM_sem = avg_animals(NREMbefore_REM_tot)
REMbefore_wake_tot, REMbefore_NREM_tot, REMbefore_REM_tot,REMbefore_totaltrials = perc_per_animal_single_stage_variable(justData, trialOccurances,'REM')
REMbefore_wake_avg, REMbefore_wake_sem = avg_animals(REMbefore_wake_tot)
REMbefore_NREM_avg, REMbefore_NREM_sem = avg_animals(REMbefore_NREM_tot)
REMbefore_REM_avg, REMbefore_REM_sem = avg_animals(REMbefore_REM_tot)
wakebefore_wake_tot, wakebefore_NREM_tot, wakebefore_REM_tot, WAKEbefore_totaltrials = perc_per_animal_single_stage_variable(justData,trialOccurances, 'Wake')
wakebefore_wake_avg, wakebefore_wake_sem = avg_animals(wakebefore_wake_tot)
wakebefore_NREM_avg, wakebefore_NREM_sem = avg_animals(wakebefore_NREM_tot)
wakebefore_REM_avg, wakebefore_REM_sem = avg_animals(wakebefore_REM_tot)
plotPercSEM(wake_avg, wake_sem, NREM_avg, NREM_sem, REM_avg, REM_sem, sec, global_name)
plotPercSEM(NREMbefore_wake_avg, NREMbefore_wake_sem, NREMbefore_NREM_avg, NREMbefore_NREM_sem, NREMbefore_REM_avg, NREMbefore_REM_sem, sec, NREM_before)
plotPercSEM(REMbefore_wake_avg, REMbefore_wake_sem, REMbefore_NREM_avg, REMbefore_NREM_sem, REMbefore_REM_avg, REMbefore_REM_sem, sec, REM_before)
plotPercSEM(wakebefore_wake_avg, wakebefore_wake_sem, wakebefore_NREM_avg, wakebefore_NREM_sem, wakebefore_REM_avg, wakebefore_REM_sem, sec, wake_before)


# trials4bootstrap = trials_by_animal(justData,trialOccurances)
# trials4bootstrap_realigned = undueBootstrapFormatting(trials4bootstrap)
# bootstrap_w, bootstrap_n, bootstrap_r = state_probabilities_all_bootstrap_better(trials4bootstrap_realigned, trialOccurances)
# wake_CI = extractConfidenceIntervals_all(bootstrap_w)
# NREM_CI = extractConfidenceIntervals_all(bootstrap_n)
# REM_CI = extractConfidenceIntervals_all(bootstrap_r)
# wake_CI = generatePlotableCI(wake_avg, wake_CI)
# NREM_CI = generatePlotableCI(NREM_avg, NREM_CI)
# REM_CI = generatePlotableCI(REM_avg, REM_CI)
# plotPercCI(wake_avg, wake_CI, NREM_avg, NREM_CI, REM_avg, REM_CI, sec, global_name +'CI')
# wake_diff, NREM_diff, REM_diff = mean_probabilites(wake_avg, NREM_avg, REM_avg)
# wake_diff_bs, NREM_diff_bs, REM_diff_bs = mean_prob_bootstrap(bootstrap_w, bootstrap_n,bootstrap_r)
# wake_inc, wake_dec = increaseORdecrease(wake_diff_bs)
# NREM_inc, NREM_dec = increaseORdecrease(NREM_diff_bs)
# REM_inc, REM_dec = increaseORdecrease(REM_diff_bs)