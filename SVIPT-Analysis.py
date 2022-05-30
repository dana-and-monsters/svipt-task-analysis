# This is the analysis script for "Transcranial direct current stimulation leads to faster acquisition of motor skills, but effects are not maintained at retention."
# A pinch force transducer was used to train participants on a motor learning task called the Sequential Visual Isometric Pinch Task (SVIPT).
# This script goes through each of the peaks in force and checks if each peak is in any window.
# Then it checks the order of the peaks to make sure they reach each one in order.

# VARIABLES TO EDIT
rootdir = '' # ** Update with your personal data directory **
maxVelocityPercent = 0.05
cutoffTrialStart = 75 # average of the first 10 trials
beginIndStart = 0
endIndStart = 10
#target windows (referencing position on the screen)
home = 50
target1start = 469.76
target1end = 541.76
target2start = 1153.4
target2end = 1225.4
target3start = 697.64
target3end = 769.64
target4start = 241.88
target4end = 313.88
# target 5 does not matter according to NK
desiredTrial = 71
# peak detection variables
threshold = 230 # 190 is below the beginning of target 4
min_distance = 50 # minimal number of samples to distinguish peak - will get smaller as subject accelerates

# Init SPSS file columns
ParticipantNum = [];SessionNum= []; BlockNum= []; Skill= []; AvgMoveTime= []; Accuracy= []; ErrorRate= [];

#### PLOT TRIALS - use this commented code to visualize trials.
            # pyplot.figure(figsize=(10,6))
            # plot(x, y, indexes) # if presmooth: write presmooth_ in front of each var
            # pyplot.title('First estimate')

## Import libraries
import numpy
import dill
import pandas
import peakutils
import os
import scipy
from scipy.signal import butter, lfilter, freqz
import smooth

from peakutils.plot import plot
from matplotlib import pyplot

class Target:
    def __init__(self):
        self.Score = []
        self.PeakValue = []
        self.PeakSample = []


for subdir, dirs, files in os.walk(rootdir):
    for file in files:
        filepath = subdir + os.sep + file
        if "warmup" in filepath: # skip warmup trials
            continue
        elif filepath.endswith(".dat"):
            with open(filepath,'r') as f:
                data = pandas.DataFrame(l.rstrip().split() for l in f)
        else:
            continue

        # initiate cols of summary data
        trialNum = []
        trialTime = []
        trialScore = []
        trialType = []
        participants = []
        days = []
        blockNum = []
        avgMoveTime = [] #list of block avg movement times (length depends on number of blocks. acquisition should = 6 blocks)
        corrTrials = []
        incorrTrials = []
        errorRate = []
        accuracy = []
        skill = []
        participantNum = []
        sessionNum = []

        #init vars
        Target1 = Target()
        Target2 = Target()
        Target3 = Target()
        Target4 = Target()
        excludedTrials = []

        subna = filepath[-11:-9]
        session = filepath[-5]
        pixeldata = data[1::4] # Select rows with data. (Every 4 rows starting at the 2nd row or ind = 1). These rows contain data representing the horizontal position of the cursor on the screen.
        pixels = numpy.array(pixeldata)
        pixeldat = pixels.astype(numpy.float)

        endofscreen = numpy.nanmax(pixeldat)
        # specify parameters for peak indexes
        thresh = threshold / endofscreen  # 190 is below the beginning of target 4

        # Analyze each trial
        trialN = 0
        for trial in pixeldat: #for row in pixeldat
            trialN += 1
            trialNum.append(trialN) # col of 1:number of trials
            # preprocessing
            presmooth_x = numpy.linspace(0, len(trial), len(trial))
            x = numpy.linspace(0, len(trial), len(trial)+10)
            presmooth_y = trial
            whereNan = numpy.isnan(presmooth_y)
            presmooth_y[whereNan] = 0
            # smooth data
            y = smooth.smooth(presmooth_y)
            presmooth_indexes = peakutils.peak.indexes(presmooth_y, thres = thresh, min_dist = min_distance) # thres = normalized threshold.
            indexes = peakutils.peak.indexes(y, thres = thresh, min_dist = min_distance) # thres = normalized threshold.
            #### PLOT TRIALS - uncomment section below to plot the trials
            # pyplot.figure(figsize=(10,6))
            # plot(x, y, indexes)
            # pyplot.title('First estimate')
            # if trialN==desiredTrial:
            #     pyplot.figure(figsize=(10, 6))
            #     plot(x, y, indexes) # if presmooth: write presmooth_ in front of each var
            target1Peaks =[]
            target1Ind = []
            target2Peaks = []
            target2Ind = []
            target3Peaks = []
            target3Ind = []
            target4Peaks = []
            target4Ind = []
            sampleList = []
            for ind in indexes:
            #check to see if any of the peaks are in any of the targets
                peak = y[ind] #pixel location
                if peak>= target1start and peak <= target1end:
                    target1Peaks.append(peak) #location along screen
                    target1Ind.append(ind) #sample at which peak took place
                elif peak>= target2start and peak <= target2end:
                    target2Peaks.append(peak)
                    target2Ind.append(ind)
                elif peak>= target3start and peak <= target3end:
                    target3Peaks.append(peak)
                    target3Ind.append(ind)
                elif peak>= target4start and peak <= target4end:
                    target4Peaks.append(peak)
                    target4Ind.append(ind)

            # Choose the maximum peak in the target as the peak (furthest point across the screen prior to cursor returning to home)
            if len(target1Peaks) > 0: # if there is at least
                # one peak in the target then its correct
                Target1.Score.append(1)
                peakIn1 = max(target1Peaks)
                peakSample1 = target1Ind[target1Peaks.index(peakIn1)]
                Target1.PeakValue.append(peakIn1) # location along screen
                Target1.PeakSample.append(peakSample1) #when target was reached
            else: # if there are no peaks, the trial is incorrect
                Target1.Score.append(0) #find peak closest to window
                Target1.PeakValue.append(0)
                Target1.PeakSample.append(0)
                peakSample1 = 0
            if len(target2Peaks) > 0:
                Target2.Score.append(1)
                peakIn2 = max(target2Peaks)
                peakSample2 = target2Ind[target2Peaks.index(peakIn2)]
                Target2.PeakValue.append(peakIn2)
                Target2.PeakSample.append(peakSample2)
            else:
                Target2.Score.append(0)
                Target2.PeakValue.append(0)
                Target2.PeakSample.append(0)
                peakSample2 = 0
            if len(target3Peaks) > 0:
                Target3.Score.append(1)
                peakIn3 = max(target3Peaks)
                peakSample3 = target3Ind[target3Peaks.index(peakIn3)]
                Target3.PeakValue.append(peakIn3)
                Target3.PeakSample.append(peakSample3)
            else:
                Target3.Score.append(0)
                Target3.PeakValue.append(0)
                Target3.PeakSample.append(0)
                peakSample3 = 0
            if len(target4Peaks) > 0:
                Target4.Score.append(1)
                peakIn4 = max(target4Peaks)
                peakSample4 = target4Ind[target4Peaks.index(peakIn4)]
                Target4.PeakValue.append(peakIn4)
                Target4.PeakSample.append(peakSample4)
            else:
                Target4.Score.append(0)
                Target4.PeakValue.append(0)
                Target4.PeakSample.append(0)
                peakSample4 = 0

            # EXCLUDE TRIALS
            # for each trial, check the order of the peaks to make sure that they reached the targets in the correct order: 1,2,3,4
            peakSampleList = [peakSample1, peakSample2, peakSample3, peakSample4] # the order in time should be 1,2,3,4
            sortedPeakSampleList = sorted(peakSampleList)
            # Exclude trial if first point is greater than 469.76
            if trial[0] > 469.76:
                trialScore.append(numpy.nan)
                excludedTrials.append(trialN)
                trialType.append('Extraordinary>469.76')

            ## April 30, 2018 - Nirsan wants to remove the criterion of checking pts go to home (start position) between trials. Uncomment below to reintroduce this criterion.
            ## Exclude trials that do not have 5 consecutive 50's of average is greater than cutoff. This was the definition for returning to the start position.
            # if smooth.check4five50s(trial, beginIndStart, endIndStart)!=5 or beginAverage > cutoffTrialStart:
            #     trialScore.append(0)
            #     excludedTrials.append(trialN)
            #     trialType.append('Excluded_Beginning')
            # elif smooth.check4five50s(trial, peakSample1, peakSample2) != 5 or smooth.check4five50s(trial, peakSample2, peakSample3) != 5 or smooth.check4five50s(trial, peakSample3, peakSample4) != 5 or smooth.check4five50s(trial, peakSample4, len(trial)-1) != 5:
            #     trialScore.append(0)
            #     excludedTrials.append(trialN)
            #     trialType.append('Incorrect_MissingReturnHome')
            # if target list in order, and there is an index in each trial, count as correct

            ## SCORING
            if sortedPeakSampleList == [peakSample1, peakSample2, peakSample3, peakSample4] and peakSample1!=0:
                trialScore.append(1)
                trialType.append('Correct')
            else:
                trialScore.append(0)
                trialType.append('IncorrectOrder')

            #### Movement time ####
            if trialType[-1]=='Extraordinary>469.76':
                trialTime.append(numpy.nan)
            else:
                firstRealPeak = indexes[numpy.argmax(y[indexes]>300)]
                if firstRealPeak <2:
                    firstRealPeak = indexes[1]
                zero2firstpeakX = range(0,firstRealPeak)
                zero2firstpeakY = y[zero2firstpeakX]
                velocity = numpy.diff(zero2firstpeakY)/numpy.diff(zero2firstpeakX)
                # peak velocity
                maxVelocity = max(velocity)

                # EDIT March 9th. NK & JC want to update movement start time as 5% of max velocity
                # Recall that maxVelocityPercent = 0.05. 5% of peak velocity = theshold for response time
                responseThresh = maxVelocity*maxVelocityPercent
                # start of movement: force_generation_onset_time = 5%peakvelocitytime
                force_generation_onset_time = numpy.argmax(velocity>responseThresh) #finds first index where true
                move_start = force_generation_onset_time

                ## EDIT January 8th: Nirsan says movement end should be the full length of trial. Uncomment the code below to define movement end as when velocity drops below 5% of max velocity.
                ## find where velocity falls below 5% of max velocity
                #flippedVelocity = numpy.flip(velocity, 0)
                #indFlippedVel= numpy.argmax(flippedVelocity>=responseThresh) #index where velocity just above threshold
                #move_end = len(velocity-indFlippedVel) #index in velocity where movement ends

                trial_length = numpy.count_nonzero(y)-15 #15 samples at end of trial
                #convert to time (sampling freq = 100Hz therefore to convert to seconds divide by 100)
                movement_time_seconds = (trial_length-move_start)/100
                trialTime.append(movement_time_seconds)

        detailed_data = numpy.array([trialNum,trialTime, trialScore, trialType, Target1.Score, Target1.PeakValue, Target1.PeakSample, Target2.Score, Target2.PeakValue, Target2.PeakSample, Target3.Score, Target3.PeakValue, Target3.PeakSample, Target4.Score, Target4.PeakValue, Target4.PeakSample]) # Error.Score, Error.PeakValue, Error.PeakSample
        dataT = detailed_data.T
        # Save matrix
        fname = 'Processed-{}.csv'.format(file[:-4])
        numpy.savetxt(fname, dataT, delimiter=",", fmt = '%s', header = 'trialNum, trialTime, trialScore, trialType, Target1.Score, Target1.PeakValue, Target1.PeakSample, Target2.Score, Target2.PeakValue, Target2.PeakSample, Target3.Score, Target3.PeakValue, Target3.PeakSample, Target4.Score, Target4.PeakValue, Target4.PeakSample', comments = '')

        totalScore = (numpy.nansum(trialScore,0))/len(trialNum)

        # ######## CALCULATE SAF ###########
        # exclude excluded trials
        class Block:
            def __init__(self,start, end):
                self.start = start
                self.end = end

        block1 = Block(1,41) # 1-40
        block2 = Block(41,71) # 41-70
        block3 = Block(71, 101) # 71-100
        block4 = Block(101,131) # 101-130
        block5 = Block(131, 161)
        block6 = Block(161,201)

        if len(pixeldat) == 40:
            AcqBlocks = [block1]
        else:
            AcqBlocks = [block1, block2, block3, block4, block5, block6]
        nBlock = 0
        for block in AcqBlocks:
            nBlock+=1
            block.averageMoveTime = numpy.nanmean(trialTime[block.start-1:block.end-1])
            b = 5.424
            block.numCorrectTrials = numpy.nansum(trialScore[block.start-1:block.end-1],0)
            c = 0
            for t in excludedTrials:
                if t in range(block.start,block.end):
                    c+=1
            block.totalTrials = (block.end-block.start)-c #c: count of excluded trials in block
            block.numIncorrectTrials = block.totalTrials-block.numCorrectTrials #total-correct = incorrect
            block.error_rate = block.numIncorrectTrials/block.totalTrials
            block.accuracy = block.numCorrectTrials/block.totalTrials
            block.skill = (1-block.error_rate)/(block.error_rate*numpy.log(block.averageMoveTime)**b)
            #make list of values for data output
            avgMoveTime.append(block.averageMoveTime)
            corrTrials.append(block.numCorrectTrials)
            incorrTrials.append(block.numIncorrectTrials)
            errorRate.append(block.error_rate)
            accuracy.append(block.accuracy)
            skill.append(block.skill)
            blockNum.append(nBlock)
            participantNum.append(subna)
            sessionNum.append(session)

            # Append to spss input vars
            ParticipantNum.append(subna)
            SessionNum.append(session)
            BlockNum.append(nBlock)
            Skill.append(block.skill)
            AvgMoveTime.append(block.averageMoveTime)
            Accuracy.append(block.accuracy)
            ErrorRate.append(block.error_rate)

        summary_data = numpy.array([blockNum,avgMoveTime,corrTrials, incorrTrials,errorRate, accuracy, skill])
        summarydataT = summary_data.T

        # Save matrix
        fname = 'SummaryData-Processed-{}.csv'.format(file[:-4])
        numpy.savetxt(fname, summarydataT, delimiter=",", fmt = '%s', header = 'Block, AverageMovementTime, CorrectTrials,IncorrectTrials, ErrorRate, Accuracy, Skill', comments = '')

spss_input_data = numpy.array([ParticipantNum,SessionNum, BlockNum, Skill, AvgMoveTime, Accuracy, ErrorRate])
spss_input_dataT = spss_input_data.T
# Save matrix
fname = 'SPSS-Input.csv'
numpy.savetxt(fname, spss_input_dataT, delimiter=",", fmt = '%s', header = 'Participant, Session, Block, Skill, AvgMovementTime, Accuracy, ErrorRate', comments = '')

