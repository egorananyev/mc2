####!/usr/bin/arch -i386 /usr/bin/python # -*- coding: utf-8 -*-
"""
Motion Clouds: SF Bandwidth (B_sf)
2016-07-29
"""

from __future__ import division  # so that 1/3=0.333 instead of 1/3=0
from psychopy import visual, core, data, event, gui
from psychopy.constants import *  # things like STARTED, FINISHED
import numpy as np
import pandas as pd
from datetime import datetime
import os, shutil, itertools, copy  # handy system and path functions
#import pyglet
#import MotionClouds as mc

#Initiating the keyboard
from psychopy.iohub import launchHubServer
io = launchHubServer()
kb_device = io.devices.keyboard

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))

# ====================================================================================
## Initial variables.
et = 0
expName = 'mc2_tgT-mcBv'
# Window circles (specified in degrees of visual angles [dva]):
#winSz = 7.2 # 5.03; calculated as 5/x=sqrt(2)/2 => x=10/sqrt(2)
winOffX = 4.25 # 6 # 5.62
winOffY = 3.5 # 5.5 (3.5cm ~= 124px)
winThickness = 2 # in pixels
# Timing variables:
ISIduration = .25
fixSz = .15
# MCs:
precompileMode = 1 # get the precompiled MCs
grtSize = 256 # size of 256 is 71mm, or 7.2dova
# Dimensions:
###### 7.2dova = 71mm = 256px; 475x296mm, 563mm viewing dist ######
dr = (1680,1050) # display resolution in px
#dr = (1366,768)
#dd = (47.5,29.6) # display dimensions in cm
dd = (29.5,16.6)
ds = 50+2.5+3.5 #49.5 # distance to screen in cm
trialNfb = False # do we give the trial number feedback?

# ====================================================================================
# Converter functions:
def cm2px(cm,dr=dr,dd=dd):
    px = int(cm*(dr[0]/dd[0]))
    return px
def px2cm(px,dr=dr,dd=dd):
    cm = px/(dr[0]/dd[0])
    return cm
def cm2dg(cm,ds=ds):
    dg = np.degrees(np.arctan(cm/ds))
    return dg
def dg2cm(dg,ds=ds):
    cm = ds*np.tan(np.radians(dg))
    return cm
def px2dg(px,cm2dg=cm2dg,px2cm=px2cm):
    dg = cm2dg(px2cm(px))
    return dg
def dg2px(dg,cm2px=cm2px,dg2cm=dg2cm):
    px = int(cm2px(dg2cm(dg)))
    return px

# ====================================================================================
# Converting win dimensions to pixels
#winSz = dg2px(winSz)
winSz = grtSize + 2
winOffX = dg2px(winOffX)
winOffY = dg2px(winOffY)
fixSz = dg2px(fixSz)
posCentL = [-winOffX, winOffY]
posCentR = [winOffX, winOffY]
#print winSz 
#print posCentL 
#print posCentR 

# ====================================================================================
# Eye tracking initialization

if et:
    import pylink as pl
    #cp = (0.4,0.4) # calibration proportion
    cd = 32

    eyeLink = ("100.1.1.1")

    displayInfo = pl.getDisplayInformation()
    print displayInfo.width, displayInfo.height
    screenCenter = (int(dr[0]/2), int(dr[1]/2))
    calScreenCenter = (int(screenCenter[0]+winOffX),
                    int(screenCenter[1]-winOffY))
    calTargDist = int(winSz/3)
    calTarg1 = calScreenCenter
    calTarg2 = (int(calScreenCenter[0]-calTargDist), int(calScreenCenter[1]))
    calTarg3 = (int(calScreenCenter[0]+calTargDist), int(calScreenCenter[1]))

    def elEndRec(el):
        # Ends the recording; adds 100ms to catch final events
        pl.endRealTimeMode()
        pl.pumpDelay(100)
        el.stopRecording()

    def eyeTrkInit (dr):
        el = pl.EyeLink()
        # sending the screen dimensions to the eye tracker:
        el.sendCommand('screen_pixel_coords = 0 0 %d %d' %dr)
        el.sendMessage('DISPLAY_COORDS 0 0 %d %d' %dr)
        el.sendCommand('generate_default_targets = NO')
        el.sendCommand('calibration_targets = %d,%d %d,%d %d,%d' % (
                       calTarg1[0], calTarg1[1],
                       calTarg2[0], calTarg2[1],
                       calTarg3[0], calTarg3[1]) )
        el.sendCommand('validation_targets = %d,%d %d,%d %d,%d' % (
                       calTarg1[0], calTarg1[1],
                       calTarg2[0], calTarg2[1],
                       calTarg3[0], calTarg3[1]) )
        # parser configuration 1 corresponds to high sensitivity to saccades:
        el.sendCommand('select_parser_configuration 1')
        # turns off "scenelink camera stuff", i.e., doesn't record the ET video
        el.sendCommand('scene_camera_gazemap = NO')
        # converting pupil area to diameter
        el.sendCommand('pupil_size_diameter = %s'%('YES'))
        return(el)
    el = eyeTrkInit(dr)
    print 'Finished initializing the eye tracker.'

    def eyeTrkCalib (el=el,dr=dr,cd=cd):
        # "opens the graphics if the display mode is not set"
        pl.openGraphics(dr,cd)
        pl.setCalibrationColors((255,255,255),(0,177,177))
        pl.setTargetSize(10, 5) 
        pl.setCalibrationSounds("","","")
        el.setCalibrationType('H3')
        pl.setDriftCorrectSounds("","off","off")
        el.disableAutoCalibration()
        el.doTrackerSetup()
        el.drawCalTarget(calTarg1)
        el.drawCalTarget(calTarg2)
        el.drawCalTarget(calTarg3)
        pl.closeGraphics()
        el.setOfflineMode()

# ====================================================================================
# Store info about the experiment session
expInfo = {u'session': u'', u'participant': u''}
dlg = gui.DlgFromDict(dictionary=expInfo, title=expName) # dialogue box
if dlg.OK == False: core.quit()  # user pressed cancel
timeNow = datetime.now()
expInfo['time'] = datetime.now().strftime('%Y-%m-%d_%H%M')
expInfo['expName'] = expName

# Setup the Window
win = visual.Window(size=dr, fullscr=True, screen=0, allowGUI=False, 
      allowStencil=False, color='grey', blendMode='avg', useFBO=True, units='pix')
# store frame rate of monitor if we can measure it successfully:
frameRate=win.getActualFrameRate()
if frameRate!=None:
    frameDur = 1.0/round(frameRate)
else:
    frameDur = 1.0/60.0 # couldn't get a reliable measure so guess

# ====================================================================================
# Eye-tracking setup

if et:
    def etSetup(el=el,dr=dr,cd=cd):
        blockLabel=visual.TextStim(win, text="Press spacebar", pos=posCentR,
                                   color="white", bold=True, alignHoriz="center",
                                   height=0.5)
        notdone=True
        while notdone:
            blockLabel.draw()
            win.flip()
            keySpace = event.getKeys(keyList=['escape','space'])
            if 'space' in keySpace:
                print 'spacebar pressed'
                eyeTrkCalib()
                win.winHandle.activate()
                print '///Finished calibration///'
                notdone=False
            elif 'escape' in keySpace:
                print 'procedure terminated'
                notdone=False
    etSetup()

    def drCor(el=el,dr=dr,cd=cd):
        pl.openGraphics(dr,cd)
        el.doDriftCorrect(calScreenCenter[0], calScreenCenter[1], 1, 0)
        pl.closeGraphics()
        print '///Finished drift correction///'

# ====================================================================================

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
if precompileMode:
    precompiledDir = '..' + os.sep + 'precompiledMCs'
dataDir = '..' + os.sep + 'data'
fileName = '%s_p%s_s%s_%s' %(expName, expInfo['participant'], expInfo['session'],
    expInfo['time'])
filePath = dataDir + os.sep + fileName
print filePath

if et:
    edfFileName = 'data.edf'
    el.openDataFile(edfFileName)
    el.sendCommand("file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,\
                    MESSAGE,BUTTON,INPUT")
    el.sendCommand("file_sample_data  = LEFT,RIGHT,GAZE,AREA,GAZERES,STATUS,\
                    HTARGET,INPUT")
    print '///set up the EDF file for eye-tracking///'

# Condition-related variables
conditionsFilePath = 'cond-files'+os.sep+'cond-'+expName+'.csv' #TEMP
print conditionsFilePath
os.chdir(_thisDir)

# ====================================================================================

endExpNow = False  # flag for 'escape' or other condition => quit the exp

# Initialize components for Routine "instructions"
instructionsClock = core.Clock()
instrTextL = visual.TextStim(win, text='Press any key to start', font='Cambria',
                             pos=posCentL, height=dg2px(.65), wrapWidth=dg2px(5),
                             color='white', alignHoriz='center')
instrTextR = visual.TextStim(win, text='Press any key to start', font='Cambria',
                             pos=posCentR, height=dg2px(.65), wrapWidth=dg2px(5),
                             color='white', alignHoriz='center')

# Initialize components for Routine "trial"
trialClock = core.Clock()
moveClock = core.Clock()
maskMoveClock = core.Clock()
ISI = core.StaticPeriod(win=win, screenHz=frameRate, name='ISI')
# circular windows:
winL = visual.Polygon(win, edges=36, size=[winSz, winSz], pos=posCentL,
                      lineWidth=winThickness, lineColor='white')
winR = visual.Polygon(win, edges=36, size=[winSz, winSz], pos=posCentR,
                      lineWidth=winThickness, lineColor='white')
# target gabor
targGab = visual.GratingStim(win, tex='sin', mask='circle', size=[winSz, winSz],
                      pos=posCentL)
# fixation:
fixL = visual.ShapeStim(win, pos=posCentL, vertices=((0,-fixSz), (0,fixSz), (0,0), 
                                                     (-fixSz,0), (fixSz,0)),
                        lineWidth=.2, closeShape=False, lineColor='white')
fixR = visual.ShapeStim(win, pos=posCentR, vertices=((0,-fixSz), (0,fixSz), (0,0), 
                                                     (-fixSz,0), (fixSz,0)),
                        lineWidth=.2, closeShape=False, lineColor='white')
# Trial number feedback:
if trialNfb:
    trialNfbText = visual.TextStim(win=win, text='', font='Cambria', 
                                   pos=(0,0), height=dg2px(.55), wrapWidth=dg2px(4.5),
                                   color='white')
# pause text:
pauseTextL = visual.TextStim(win, text='Press Spacebar to continue', font='Cambria',
                             alignHoriz='center', pos=posCentL, height=dg2px(.7),
                             wrapWidth=dg2px(3), color='white')
pauseTextR = visual.TextStim(win, text='Press Spacebar to continue', font='Cambria',
                             alignHoriz='center', pos=posCentR, height=dg2px(.7),
                             wrapWidth=dg2px(3), color='white')

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 

#------Prepare to start Routine "instructions"-------
t = 0
instructionsClock.reset()  # clock 
frameN = -1
# update component parameters for each repeat
instrKey = event.BuilderKeyResponse()  # create an object of type KeyResponse
instrKey.status = NOT_STARTED
# keep track of which components have finished
instructionsComponents = []
instructionsComponents.append(instrTextL)
instructionsComponents.append(instrTextR)
instructionsComponents.append(instrKey)
for thisComponent in instructionsComponents:
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED


# ====================================================================================
# Setting up the conditions:
condList = data.importConditions(conditionsFilePath)
stairs = []
for thisCond in condList:
    thisInfo = copy.copy(thisCond)
    stairLabel = 'st-' + str(thisCond['startContr']) + \
                 '_mcBv-' + str(thisCond['mcBv']) + \
                 '_targTpeak-' + str(thisCond['targTpeak'])
    thisInfo['label'] = stairLabel
    nTrials = thisCond['trialN']
    thisStair = data.QuestHandler(startVal = thisInfo['startContr'],
                                  extraInfo = thisInfo,
                                  startValSd = .2, pThreshold = .63,
                                  gamma = 0.01, nTrials = nTrials,
                                  minVal=0, maxVal=1)
    stairs.append(thisStair)

# An empty data set for storing behavioural responses:
behResp = []
    
# Creating a copy of the Conditions file for book-keeping and analyses:
if not os.path.exists(filePath):
    os.makedirs(filePath)
shutil.copyfile(conditionsFilePath, filePath + os.sep + 
                os.path.basename(conditionsFilePath))
dataFileName = filePath + os.sep + fileName + '.csv'

# ====================================================================================
# Various functions for use in trials:

def sigmoid(x):
    return 1 / (1 + np.exp(-x))

def sigmoidMod(x): # modified such that 0->0, 1->1
    return 1 / (1 + np.exp(-x*10+5))

x = np.arange(-grtSize/2,grtSize/2)
y = np.arange(-grtSize/2,grtSize/2)
x, y = np.meshgrid(x, y)
R = np.sqrt((x+.5)**2 + (y+.5)**2) # adding .5 ensures symmetry

def periMask(periGap, periFade, R=R):
    return sigmoid(R * (-10./(periFade)) + 5 + periGap*(10./periFade))*2 - 1

# ====================================================================================

#-------Start Routine "instructions"-------
continueRoutine = True
while continueRoutine:
    # get current time
    t = instructionsClock.getTime()
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *instrText* updates
    if t >= 0.0 and instrTextL.status == NOT_STARTED:
        # keep track of start time/frame for later
        instrTextL.tStart = t  # underestimates by a little under one frame
        instrTextL.frameNStart = frameN  # exact frame index
        instrTextL.setAutoDraw(True)
        instrTextR.tStart = t  # underestimates by a little under one frame
        instrTextR.frameNStart = frameN  # exact frame index
        instrTextR.setAutoDraw(True)
    
    # *instrKey* updates
    if t >= 0.0 and instrKey.status == NOT_STARTED:
        # keep track of start time/frame for later
        instrKey.tStart = t  # underestimates by a little under one frame
        instrKey.frameNStart = frameN  # exact frame index
        instrKey.status = STARTED
        # keyboard checking is just starting
        event.clearEvents(eventType='keyboard')
        winL.setAutoDraw(True)
        winR.setAutoDraw(True)
    if instrKey.status == STARTED:
        theseKeys = event.getKeys()
        
        # check for quit:
        if "escape" in theseKeys:
            endExpNow = True
        if len(theseKeys) > 0:  # at least one key was pressed
            # a response ends the routine
            continueRoutine = False
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        routineTimer.reset()  # if we abort early the non-slip timer needs reset
        break
    continueRoutine = False  # reverts to True if at least 1 component still running
    for thisComponent in instructionsComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # check for quit (the Esc key)
    if endExpNow or event.getKeys(keyList=["escape"]):
        core.quit()
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()
    else:  # this Routine was not non-slip safe so reset non-slip timer
        routineTimer.reset()

#-------Ending Routine "instructions"-------
for thisComponent in instructionsComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)

# ====================================================================================
# Initiating the trial loop

nDone=0
for trialN in range(nTrials):
    np.random.shuffle(stairs)
    for thisStair in stairs:

        nDone += 1
        if trialNfb:
            trialNfbText.text = str(nDone) + '/' + str(trials.nTotal)
        trialNstr = '#' + str(nDone) + '/' + str(nTrials*len(stairs))

        ## Setting up trial variables

        # current contrast:
        thisContr = thisStair.next()
        contrStr = 'start=%.1f, cur=%.2f' %(thisStair.extraInfo['startContr'], thisContr)

        # mc mask:
        mcSz = thisStair.extraInfo['mcSz']
        mcSf = thisStair.extraInfo['mcSf']
        mcBv = thisStair.extraInfo['mcBv']
        mcBsf = thisStair.extraInfo['mcBsf']

        # target:
        targSz = thisStair.extraInfo['targSz']
        targSf = thisStair.extraInfo['targSf']
        targXoff = thisStair.extraInfo['targXoff']
        targYoff = thisStair.extraInfo['targYoff']
        targV = thisStair.extraInfo['targV']

        # Target response criteria:
        targOri1 = thisStair.extraInfo['targOri1']
        targOri2 = thisStair.extraInfo['targOri2']
        if not targOri1 == targOri2: # if the two targ oris are not the same, decided randomly
            allTargOris = np.array([targOri1, targOri2])
            thisTargOri = allTargOris[np.random.randint(2)]

        # Setting up the target with the above characteristics:
        targGab.size = targSz
        targGab.sf = targSf
        targGab.pos = targGab.pos + np.array([targXoff, targYoff])
        targGab.ori = thisTargOri

        # Temporal variables:
        targTtot = thisStair.extraInfo['targTtot']
        targTpeak = thisStair.extraInfo['targTpeak']
        targTstart = targTpeak-(targTtot/2)
        targTend = targTpeak+(targTtot/2)
        trialT = thisStair.extraInfo['trialT'] # -win.monitorFramePeriod*0.75
        
        print 'TRIAL' + '\t' + 'CONTRAST' + '\t\t' + 'mcBv' + '\t' + 'targOri' + '\t' + 'targTpeak'
        print trialNstr + '\t' + contrStr + '\t' + str(mcBv) + '\t' + str(thisTargOri) + '\t' + str(targTpeak)

        # view setup: fade, gap, and fixation cross
        fixCross = thisStair.extraInfo['fixCross']
        periFade = thisStair.extraInfo['mcPeriFade']
        periGap = np.int(mcSz / 2 - periFade)
        #print 'periFade=' + str(periFade) + '; periGap=' + str(periGap)

        nFrames = 60 # number of frames per sequence
        
        # creating an empty matrix for keeping the behavioural responses:
        behRespTrial = []
            
        # initiating the mc gratings:
        mcV = 0
        grt = np.load(precompiledDir + os.sep + 'mc_' + str(mcV) +
                '_sf' + str(mcSf) + '_bsf' + str(mcBsf) + '_bv' + str(mcBv) + 
                '_sz' + str(mcSz) + '.npy')

        # creating a mask, which is fixed for a given trial:
        mcPeriMask = periMask(periGap, periFade)

        #------prepare to start routine "trial"-------
        t = 0
        trialClock.reset()  # clock 
        frameN = -1

        # anchors:
        elStopped = False
        keyPause = False
        targRespGiven = False
        behRespRecorded = False

        # update component parameters for each repeat
        key_arrow = event.BuilderKeyResponse()  # create an object of type keyresponse
        key_arrow.status = NOT_STARTED
        # keep track of which components have finished
        trialComponents = []
        trialComponents.append(winL)
        trialComponents.append(winR)
        trialComponents.append(pauseTextL)
        trialComponents.append(pauseTextR)
        trialComponents.append(fixL)
        trialComponents.append(fixR)
        trialComponents.append(key_arrow)
        trialComponents.append(ISI)
        for thisComponent in trialComponents:
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        
        # ////////////////////////////////////////////////////////////////////////////////
        if et:
            el.sendMessage("TRIALID " + str(nDone))
            trialStartStr = datetime.now().strftime('%Y-%m-%d_%H%M%S')
            el.sendMessage("TIMESTAMP " + trialStartStr)
            el.setOfflineMode()
            pl.msecDelay(50) 
            error = el.startRecording(1,1,1,1)
        # ////////////////////////////////////////////////////////////////////////////////
        
        #-------Start Routine "trial"-------
        continueRoutine = True
        while continueRoutine:
            # get current time
            t = trialClock.getTime()
            frameN = frameN + 1 # number of completed frames (0 is the first frame)
            # update/draw components on each frame
            
            # *winL* updates
            if winL.status == NOT_STARTED:
                # keep track of start time/frame for later
                winL.tStart = t  # underestimates by a little under one frame
                winL.frameNStart = frameN  # exact frame index
                winL.setAutoDraw(True)
                winL.status = STARTED
            
            # *winR* updates
            if winR.status == NOT_STARTED:
                # keep track of start time/frame for later
                winR.tStart = t  # underestimates by a little under one frame
                winR.frameNStart = frameN  # exact frame index
                winR.setAutoDraw(True)
                winR.status = STARTED

            if trialNfb:
                trialNfbText.draw()

            # mcMask and targ presentation:
            if t < trialT:
                # mcMask:
                mcMask = visual.GratingStim(win, tex=grt[:,:,frameN%nFrames], 
                    size=(grtSize,grtSize), pos=[winOffX, winOffY], 
                    interpolate=False, mask=mcPeriMask)
                mcMask.draw()
                # target presentation:
                if t > targTstart and t < targTpeak:
                    targGab.opacity = sigmoidMod((t-targTstart)*2/targTtot)*thisContr
                elif t > targTpeak and t < targTend:
                    targGab.opacity = sigmoidMod((targTend-t)*2/targTtot)*thisContr
                else:
                    targGab.opacity = 0
                targGab.draw()
                
                # drawing the fixation cross, if any:
                if fixCross:
                    fixL.draw()
                    fixR.draw()
            
            # *key_arrow* updates for target reponses:
            if key_arrow.status == NOT_STARTED:
                # keep track of start time/frame for later
                key_arrow.tStart = t  # underestimates by a little under one frame
                key_arrow.frameNStart = frameN  # exact frame index
                key_arrow.status = STARTED
                # keyboard checking is just starting
                key_arrow.clock.reset()  # now t=0
                event.clearEvents(eventType='keyboard')
                kb_device.clearEvents()
            # registering response at throughout the trial
            if key_arrow.status == STARTED:
                theseKeys = event.getKeys(keyList=['left','right'])
                if len(theseKeys) > 0:
                    if 'left' in theseKeys:
                        print 'response: left tilt'
                        behRespTrial = targOri1 # the first number is negative
                        targRespGiven = True
                    elif 'right' in theseKeys:
                        print 'response: right tilt'
                        behRespTrial = targOri2 # the second number is positive
                        targRespGiven = True
                    if targRespGiven:
                        if behRespTrial == thisTargOri: corrResp = 1
                        else: corrResp = 0

            if t > trialT and not elStopped:
                # stopping eye-tracking recording:
                if et:
                    elEndRec(el)
                    elStopped = True

            # pause text and data exporting
            if targRespGiven and not keyPause and t>trialT:
                pauseTextL.draw()
                pauseTextR.draw()
                if 'space' in event.getKeys(keyList=['space']):
                    keyPause = True
                    # only update the response upon pressing 'space':
                    if corrResp: print 'correct'
                    else: print 'incorrect'
                    thisStair.addResponse(corrResp) 
                    #print 'spacebar pressed'
                    pauseTextL.setAutoDraw(False)
                    pauseTextR.setAutoDraw(False)

            # *ISI* period
            if ISI.status == NOT_STARTED and t>=trialT and keyPause:
                # keep track of start time/frame for later
                ISI.tStart = t  # underestimates by a little under one frame
                ISI.frameNStart = frameN  # exact frame index
                fixL.setAutoDraw(True)
                fixR.setAutoDraw(True)
                ISI.start(ISIduration)
            #one frame should pass before updating params and completing
            elif ISI.status == STARTED and t >= (ISI.tStart + ISIduration): 
                fixL.setAutoDraw(False)
                fixR.setAutoDraw(False)
                ISI.complete() #finish the static period
                continueRoutine = False
            
            # check if all components have finished
            # a component has requested a forced-end of Routine:
            if not continueRoutine: 
                # if we abort early the non-slip timer needs reset:
                routineTimer.reset() 
                break
            # will revert to True if at least one component still running
            continueRoutine = False  
            for thisComponent in trialComponents:
                if hasattr(thisComponent, "status") and \
                        thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # check for quit (the Esc key)
            if endExpNow or event.getKeys(keyList=["escape"]):
                print np.shape(behResp)
                if et:
                    elEndRec(el)
                core.quit()
            
            # refresh the screen
            # don't flip if this routine is over or we'll get a blank screen
            if continueRoutine:  
                win.flip()
            else: # this Routine was not non-slip safe so reset non-slip timer
                routineTimer.reset()
        
        #-------Ending Routine "trial"-------
        for thisComponent in trialComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)

    # thisExp.nextEntry()

nStairsDone = 0
for thisStair in stairs:
    nStairsDone += 1
    stairFileName = filePath + os.sep + thisStair.extraInfo['label']
    print thisStair.mean()
    thisStair.saveAsPickle(stairFileName)
    thisStair.saveAsText(stairFileName)
    mcPeriFade = thisStair.extraInfo['mcPeriFade']
    mcPeriGap = np.int(thisStair.extraInfo['mcSz']/2-mcPeriFade)
    # have the information recorded in a csv file as well:
    dT = pd.DataFrame({'expName': expName, 'time': expInfo['time'],
                       'participant': expInfo['participant'],
                       'session': expInfo['session'],
                       'nTrials': nTrials,
                       'mcSz': thisStair.extraInfo['mcSz'],
                       'mcSf': thisStair.extraInfo['mcSf'],
                       'mcBv': thisStair.extraInfo['mcBv'],
                       'mcBsf': thisStair.extraInfo['mcBsf'],
                       'mcPeriGap': mcPeriGap,
                       'mcPeriFade': mcPeriFade,
                       'targSz': thisStair.extraInfo['targSz'],
                       'targSf': thisStair.extraInfo['targSf'],
                       'targOri1': thisStair.extraInfo['targOri1'],
                       'targOri2': thisStair.extraInfo['targOri2'],
                       'targXoff': thisStair.extraInfo['targXoff'],
                       'targYoff': thisStair.extraInfo['targYoff'],
                       'targV': thisStair.extraInfo['targV'],
                       'targTtot': thisStair.extraInfo['targTtot'],
                       'targTpeak': thisStair.extraInfo['targTpeak'],
                       'trialT': thisStair.extraInfo['trialT'],
                       'fixCross': thisStair.extraInfo['fixCross'],
                       'stairLabel': thisStair.extraInfo['label'],
                       'stairStart': thisStair.extraInfo['startContr'],
                       'stairMean': [thisStair.mean()]})
    # to preserve the column order:
    dataCols = ['expName', 'time', 'participant', 'session', 'nTrials',
                'mcSz', 'mcSf', 'mcBv', 'mcBsf', 'mcPeriGap', 'mcPeriFade', 
                'targSz', 'targSf', 'targOri1', 'targOri2', 'targXoff',
                'targYoff', 'targV', 'targTtot', 'targTpeak', 'trialT',
                'fixCross', 'stairLabel', 'stairStart', 'stairMean']
    if nStairsDone == 1: df = dT
    else: df = pd.concat([df,dT])
# Recording the data to a csv file:
df.to_csv(dataFileName, index=False, columns=dataCols)
print 'wrote the data set to ' + dataFileName

if et:
    # File transfer and cleanup!
    pl.endRealTimeMode()
    el.setOfflineMode()						  
    pl.msecDelay(600) 

    #Close the file and transfer it to Display PC
    el.closeDataFile()
    el.receiveDataFile(edfFileName, edfFileName)
    os.rename(edfFileName, filePath + os.sep + edfFileName)
    el.close()

print "finished the experiment"

win.close()
core.quit()
