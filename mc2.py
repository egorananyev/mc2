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
import os, shutil, itertools  # handy system and path functions
#import pyglet
import MotionClouds as mc

#Initiating the keyboard
from psychopy.iohub import launchHubServer
io = launchHubServer()
kb_device = io.devices.keyboard

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))

# ====================================================================================
## Initial variables.
et = 0
expName = 'mc2_targOnXmaskBv'
# Window circles (specified in degrees of visual angles [dva]):
#winSz = 7.2 # 5.03; calculated as 5/x=sqrt(2)/2 => x=10/sqrt(2)
winOffX = 4.25 # 6 # 5.62
winOffY = 3.5 # 5.5 (3.5cm ~= 124px)
winThickness = 2 # in pixels
fdbkLen = .25 # the length of the feedback line, in degrees
fdbkThick = 4 # the tickness of the feedback line, in pixels
# Timing variables:
ISIduration = .5
fixSz = .15
# MCs:
precompileMode = 1 # get the precompiled MCs
grtSize = 256 # size of 256 is 71mm, or 7.2dova
defAlpha = .2
# Ring steps (for ct):
ringSteps = 10
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
fdbkLen = dg2px(fdbkLen)
fixSz = dg2px(fixSz)
posCentL = [-winOffX, winOffY]
posCentR = [winOffX, winOffY]
print winSz 
print posCentL 
print posCentR 
print fdbkLen 

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
expInfo = {u'session': u'', u'participant': u'', u'sat': 0}
dlg = gui.DlgFromDict(dictionary=expInfo, title=expName) # dialogue box
if dlg.OK == False: core.quit()  # user pressed cancel
timeNow = datetime.now()
expInfo['time'] = datetime.now().strftime('%Y-%m-%d_%H%M')
expInfo['expName'] = expName
sat = expInfo['sat'] # saturation for color 'red'

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
# color masks:
colMaskL = visual.GratingStim(win, size=[grtSize, grtSize], pos=posCentL, opacity=defAlpha,
                              colorSpace='hsv')
colMaskR = visual.GratingStim(win, size=[grtSize, grtSize], pos=posCentR, opacity=defAlpha,
                              colorSpace='hsv')
# direction feedback:
dirFdbkL = visual.Line(win, start=[0,0], end=[0,0], lineColor='white',
                        lineWidth=fdbkThick)
dirFdbkR = visual.Line(win, start=[0,0], end=[0,0], lineColor='white',
                        lineWidth=fdbkThick)
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
# question text:
qntxtL = visual.TextStim(win=win,
                         text='1=not stable\n2=not very stable\n3=almost stable\
                               \n4=completely stable', font='Cambria', 
                         pos=posCentL, height=dg2px(.55), wrapWidth=dg2px(4.5),
                         color='white')
qntxtR = visual.TextStim(win=win,
                         text='1=not stable\n2=not very stable\n3=almost stable\
                               \n4=completely stable', font='Cambria', 
                         pos=posCentR, height=dg2px(.55), wrapWidth=dg2px(4.5),
                         color='white')
# feedback ring (for ct):
ringL = visual.Polygon(win, edges=36, size=[winSz, winSz], ori=0, 
                       pos=posCentL, lineWidth=winThickness, lineColor='red',
                       opacity=.1, interpolate=True)
ringR = visual.Polygon(win, edges=36, size=[winSz, winSz], ori=0, 
                       pos=posCentR, lineWidth=winThickness, lineColor='red',
                       opacity=.1, interpolate=True)
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
conds = []
commonNTrials = []
for thisCondition in condList:
    nTrials = thisCondition['trialN']
    # print 'Number of trials in this condition: ' + str(nTrials)
    conds.append(thisCondition)
    commonNTrials = nTrials

# An empty data set for storing behavioural responses:
behResp = []
    
# Printing the attributes of the conds:  
print commonNTrials
trials = data.TrialHandler(conds, commonNTrials, extraInfo=expInfo)
# Creating a copy of the Conditions file for book-keeping and analyses:
if not os.path.exists(filePath):
    os.makedirs(filePath)
shutil.copyfile(conditionsFilePath, filePath + os.sep + 
                os.path.basename(conditionsFilePath))
dataFileName = filePath + os.sep + fileName + '.csv'

# ====================================================================================
# Various functions for use in trials:

# Increase/decrease of the after-trial central motion task:
def ringSzFn(ring, periGap, ringSzMulti):
    sz = ring.size[1] + ringSzMulti * (periGap*2) / ringSteps
    if sz > (periGap*2):
        sz = (periGap*2)
    elif sz < (periGap*2)/ringSteps:
        sz = (periGap*2)/ringSteps
    ring.setSize([sz,sz])
    return ring

def drawFdbkAngle(dirFdbk, lr, angle, winOffX=winOffX, winOffY=winOffY, winSz=winSz):
    lineStart = [int( lr*winOffX + np.cos(np.radians(angle))*winSz/2 ),
                 int( winOffY + np.sin(np.radians(angle))*winSz/2 )]
    lineEnd = [int( lr*winOffX + np.cos(np.radians(angle)) * (winSz/2 + fdbkLen) ),
               int( winOffY + np.sin(np.radians(angle)) * (winSz/2 + fdbkLen) )]
    dirFdbk.start = lineStart
    dirFdbk.end = lineEnd
    return dirFdbk

def sigmoid(x):
    return 1 / (1 + np.exp(-x))

x = np.arange(-grtSize/2,grtSize/2)
y = np.arange(-grtSize/2,grtSize/2)
x, y = np.meshgrid(x, y)
R = np.sqrt((x+.5)**2 + (y+.5)**2) # adding .5 ensures symmetry

def combinedMask(fovGap, fovFade, periGap, periFade, annuR, annuWidth, R=R):
    if fovFade > 0:
        fovMask = sigmoid(R * (10./fovFade) - (fovGap * (10./fovFade)) - 5)*2 - 1
    if periFade > 0:
        periMask = sigmoid(R * (-10./(periFade)) + 5 + periGap*(10./periFade))*2 - 1
    if annuR > 0:
        annuFadeOut = sigmoid(R * (-10./(annuWidth)) + 5 + annuR*(10./annuWidth))*2 - 1
        annuFadeIn = sigmoid(R * (10./annuWidth) - (annuR * (10./annuWidth)) - 5)*2 - 1
        annuMask = (annuFadeOut * annuFadeIn)*(-2)-1
    if fovFade == 0 and periFade >0 and annuR == 0:
        return periMask
    if fovFade > 0 and periFade > 0 and annuR == 0:
        return fovMask * periMask
    if fovFade == 0 and periFade > 0 and annuR > 0:
        return periMask * annuMask
    if fovFade > 0 and periFade > 0 and annuR > 0:
        return fovMask * periMask * annuMask

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

#win.winHandle.minimize()
#win.flip()
#drCor()
#win.winHandle.maximize()
#win.flip()
#win.winHandle.activate()

# ====================================================================================
# Initiating the trial loop

nDone=0
for thisTrial in trials:
    print '===new=trial==='
    nDone += 1
    if trialNfb:
        trialNfbText.text = str(nDone) + '/' + str(trials.nTotal)
    print 'trial#' + str(nDone) + '/' + str(trials.nTotal)

    ## Setting up trial variables

    # Motion setup for the trial:
    dirL = thisTrial['dirL']
    dirR = thisTrial['dirR']
    vL = thisTrial['vL']
    vR = thisTrial['vR']
    szL = thisTrial['szL']
    szR = thisTrial['szR']
    sfL = thisTrial['sfL']
    sfR = thisTrial['sfR']
    BvL = thisTrial['BvL']
    BvR = thisTrial['BvR']
    if BvL == 'NA': BvL = .5 # default value
    if BvR == 'NA': BvR = .5
    BsfL = thisTrial['BsfL']
    BsfR = thisTrial['BsfR']
    thL = thisTrial['thL']
    thR = thisTrial['thR']
    trialT = thisTrial['trialT'] # -win.monitorFramePeriod*0.75

    # Print out trial info based on the nature of the experiment:
    if not expName == 'mcEcc_ct-tXbv' and not expName == 'mcEcc_ct-szXbv':
        print 'dirL=' + str(dirL) + '; dirR=' + str(dirR)
    if expName == 'mcEcc_ct-sfXv' or expName == 'mcvct':
        print 'vL=' + str(vL) + '; vR=' + str(vR)
    if expName == 'mcEcc_ct-sfXv':
        print 'sfL=' + str(sfL) + '; sfR=' + str(sfR)
    if expName == 'mcEcc_ct-sfXbv' or expName == 'mcEcc_ct-bv':
        print 'BvL=' + str(BvL) + '; BvR=' + str(BvR)
    if expName == 'mcEcc_ct-bsfXv':
        print 'BsfL=' + str(BsfL) + '; BsfR=' + str(BsfR)
    if expName == 'mcEcc_ct-t' or expName == 'mcEcc_ct-tXbv':
        print 'trialT=' + str(trialT)
    if expName == 'mcEcc_ct-th':
        print 'thL=' + str(thL) + '; thR=' + str(thR)
    if expName == 'mcEcc_ct-offXbv':
        offX = thisTrial['offX']
        offY = thisTrial['offY']
        print 'offX=' + str(offX) + '; offY=' + str(offY)
    else:
        offX = 0; offY = 0
    if expName == 'mcEcc_ct-szRelXbv' or expName == 'mcEcc_ct-szRel0Xbv':
        szRelL = thisTrial['szRelL']
        szRelR = thisTrial['szRelR']
        print 'szRelL=' + str(szRelL) + '; szRelR=' + str(szRelR)
    else:
        szRelL = 1; szRelR = 1
    if expName == 'mcEcc_ct-tRelXbv':
        tOffL = thisTrial['tOffL']
        tOffR = thisTrial['tOffR']
        print 'tOffL=' + str(tOffL) + '; tOffR=' + str(tOffR)
    else:
        tOffL = 0; tOffR = 0

    # View setup: Fade, gap, and fixation cross
    centTask = thisTrial['centTask']
    fixCross = thisTrial['fixCross']
    fovGap = thisTrial['fovGap']
    fovFade = thisTrial['fovFade']
    periGap = thisTrial['periGap']
    periFade = thisTrial['periFade']
    #print 'fovGap=' + str(fovGap) + '; fovFade=' + str(fovFade)
    #print 'periGap=' + str(periGap) + '; periFade=' + str(periFade)
    #print 'annuR=' + str(annuR) + '; annuWidth=' + str(annuWidth)
    annuR = thisTrial['annuR']
    annuWidth = thisTrial['annuWidth']
    stabQn = thisTrial['stabQn'] # do we need to ask the qn on rivalry stability?
    # If central task, the ring size is set to twice the periGap (the actual stim size):
    if centTask:
        if expName == 'mcEcc_ct-szRelXbv':
            ringSz = int(np.max([periGap*szRelL*2,periGap*szRelR*2]))
            print 'szRelL=' + str(periGap*szRelL*2) + '; szRelR=' + str(periGap*szRelR*2)
        else:
            ringSz = int(periGap*2) # gets overwritten later if changed
        if expName == 'mcEcc_ct-szXv' or expName == 'mcEcc_ct-szXbv':
            print 'stimSz=' + str(ringSz)
        ringL.setSize([ringSz,ringSz])
        ringR.setSize([ringSz,ringSz])

    # Color, if any:
    colorEither = [[150,1,1],[330,sat,1]] # green and magenta
    if thisTrial['colDirL'] == 'rand': # picking one at random
        colorPick = np.random.permutation(colorEither)
        colorL = colorPick[0]
        colorR = colorPick[1]
    elif thisTrial['colDirL'] == thisTrial['colDirR']:
        colorL = [0,0,0] # greyscale (black)
        colorR = [0,0,0]
    else:
        # color for eye <- color for direction:
        if dirL == 180 and dirR == 0:
            colorL = colorEither[thisTrial['colDirL']] # if colDirL==0, green [0], magenta if 1
            colorR = colorEither[thisTrial['colDirR']]
        if dirL == 0 and dirR == 180:
            colorL = colorEither[thisTrial['colDirR']]
            colorR = colorEither[thisTrial['colDirL']]
    #print 'colorL = ' + str(colorL) + '; colorR = ' + str(colorR)

    nFrames = 60 # number of frames per sequence
    
    # Creating an empty matrix for keeping the behavioural responses:
    if trialT <= 3: # if short trial, response to be given at the end
        behRespTrial = []
    else: # for longer trials, enable continuous tracking:
        behRespTrial = np.empty([1, trialT*nFrames]) 
        behRespTrial[:] = np.NAN
        
    # initiating the gratings
    if precompileMode:
        grtL = np.load(precompiledDir + os.sep + 'mc_' + str(vL) +
               '_sf' + str(sfL) + '_bsf' + str(BsfL) + '_bv' + str(BvL) + 
               '_sz' + str(szL) + '.npy')
        grtR = np.load(precompiledDir + os.sep + 'mc_' + str(vR) +
               '_sf' + str(sfR) + '_bsf' + str(BsfR) + '_bv' + str(BvR) +
               '_sz' + str(szR) + '.npy')
    else:
        fx, fy, ft = mc.get_grids(szL, szL, nFrames)
        grtCol = mc.envelope_color(fx, fy, ft)
        grtL = 2*mc.rectif(mc.random_cloud(grtCol * 
               mc.envelope_gabor(fx, fy, ft, sf_0=sfL, B_sf=BsfL, B_V=BvL,
               V_X=vL, B_theta=np.inf))) - 1
        fx, fy, ft = mc.get_grids(szR, szR, nFrames)
        grtCol = mc.envelope_color(fx, fy, ft)
        grtR = 2*mc.rectif(mc.random_cloud(grtCol * 
               mc.envelope_gabor(fx, fy, ft, sf_0=sfR, B_sf=BsfR, B_V=BvR,
               V_X=vR, B_theta=np.inf))) - 1

    # Creating a mask, which is fixed for a given trial:
    if szRelL == 0:
        curMaskR = combinedMask(fovGap, fovFade, periGap*szRelR, periFade, annuR, annuWidth)
        curMaskL = curMaskR
    elif szRelR == 0:
        curMaskL = combinedMask(fovGap, fovFade, periGap*szRelL, periFade, annuR, annuWidth)
        curMaskR = curMaskL
    else:
        curMaskL = combinedMask(fovGap, fovFade, periGap*szRelL, periFade, annuR, annuWidth)
        curMaskR = combinedMask(fovGap, fovFade, periGap*szRelR, periFade, annuR, annuWidth)

    # Using the mask to assign both the greyscale values and the mask for our color masks:
    if not colorL == colorR: # same colors mean no color mask
        colMaskL.tex = (curMaskL + 1) / 2
        colMaskL.color = colorL
        colMaskL.mask = curMaskL
        colMaskR.tex = (curMaskR + 1) / 2
        colMaskR.color = colorR
        colMaskR.mask = curMaskR

    #------Prepare to start Routine "trial"-------
    t = 0
    trialClock.reset()  # clock 
    frameN = -1
    tMaskMove = 0
    qnResp = 0
    elStopped = False
    ringDrawn = False
    key_pressed = False
    key_pause = False
    respGiven = False # will be True once all questions are answered
    if not stabQn and trialT > 5: # for longer trials with continous monitoring, no need to ask the qn:
        key_qn = True
    else: # for shorter trials, or trials that have stabQn explicitely enforced, ask:
        key_qn = False
    behRespRecorded = False
    someKeyPressed = False # to prevent recording key releases at trial beginning
    # update component parameters for each repeat
    key_arrow = event.BuilderKeyResponse()  # create an object of type KeyResponse
    key_arrow.status = NOT_STARTED
    # keep track of which components have finished
    trialComponents = []
    trialComponents.append(winL)
    trialComponents.append(winR)
    trialComponents.append(dirFdbkL)
    trialComponents.append(dirFdbkR)
    trialComponents.append(qntxtL)
    trialComponents.append(qntxtR)
    if centTask:
        trialComponents.append(ringL)
        trialComponents.append(ringR)
    trialComponents.append(key_arrow)
    trialComponents.append(pauseTextL)
    trialComponents.append(pauseTextR)
    trialComponents.append(ISI)
    trialComponents.append(fixL)
    trialComponents.append(fixR)
    for thisComponent in trialComponents:
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    
    # ////////////////////////////////////////////////////////////////////////////////
    #win.winHandle.minimize()
    #drCor(el,dr,cd)
    #win.winHandle.maximize()
    #win.winHandle.activate()
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

        # stimulus presentation:
        if t < trialT:
            if t > tOffL and not szRelL == 0:
                stimL = visual.GratingStim(win, tex=grtL[:,:,frameN%nFrames], 
                    size=(grtSize,grtSize), pos=[-winOffX+dg2px(offX), winOffY+dg2px(offY)], 
                    interpolate=False, mask=curMaskL, ori=90+dirL)
                stimL.draw()
            if t > tOffR and not szRelR == 0:
                stimR = visual.GratingStim(win, tex=grtR[:,:,frameN%nFrames], 
                    size=(grtSize,grtSize), pos=[winOffX+dg2px(offX), winOffY+dg2px(offY)], 
                    interpolate=False, mask=curMaskR, ori=90+dirR)
                stimR.draw()
            # Drawing the color masks:
            if not colorL == colorR: # same colors mean no color mask
                colMaskL.draw()
                colMaskR.draw()
            # Drawing the fixation cross, if any:
            if fixCross:
                fixL.draw()
                fixR.draw()
        
        # *key_arrow* updates
        if key_arrow.status == NOT_STARTED:
            # keep track of start time/frame for later
            key_arrow.tStart = t  # underestimates by a little under one frame
            key_arrow.frameNStart = frameN  # exact frame index
            key_arrow.status = STARTED
            # keyboard checking is just starting
            key_arrow.clock.reset()  # now t=0
            event.clearEvents(eventType='keyboard')
            kb_device.clearEvents()
        # registering the response continuously:
        if key_arrow.status == STARTED and trialT > 5 and t < trialT:
            thesePresses = kb_device.getPresses(keys=['left','right','up','down'])
            theseReleases = kb_device.getReleases(keys=['left','right','up','down'])
            if len(thesePresses) > 0:
                dirFdbkL.setAutoDraw(True)
                dirFdbkR.setAutoDraw(True)
                keyPressFN = frameN
                someKeyPressed = True
                if 'left' in thesePresses:
                    print '"left" key is pressed'
                    drawFdbkAngle(dirFdbkL, -1, 180)
                    drawFdbkAngle(dirFdbkR, 1, 180)
                    whichKeyPressed = 'left' # only needed for final key press
                elif 'right' in thesePresses:
                    print '"right" key is pressed'
                    drawFdbkAngle(dirFdbkL, -1, 0)
                    drawFdbkAngle(dirFdbkR, 1, 0)
                    whichKeyPressed = 'right'
                elif 'up' in thesePresses:
                    print '"up" key is pressed'
                    drawFdbkAngle(dirFdbkL, -1, 90)
                    drawFdbkAngle(dirFdbkR, 1, 90)
                    whichKeyPressed = 'up'
                elif 'down' in thesePresses:
                    print '"down" key is pressed'
                    drawFdbkAngle(dirFdbkL, -1, 270)
                    drawFdbkAngle(dirFdbkR, 1, 270)
                    whichKeyPressed = 'down'
            if len(theseReleases) > 0 and someKeyPressed:
                dirFdbkL.setAutoDraw(False)
                dirFdbkR.setAutoDraw(False)
                print '...released after ' + \
                      str(np.around((frameN-keyPressFN)/60,2)) + 's'
                someKeyPressed = False
                if 'left' in theseReleases:
                    behRespTrial[0,keyPressFN:frameN+1] = 180 # left
                elif 'right' in theseReleases:
                    behRespTrial[0,keyPressFN:frameN+1] = 0 # right
                elif 'up' in theseReleases:
                    behRespTrial[0,keyPressFN:frameN+1] = 90 # right
                elif 'down' in theseReleases:
                    behRespTrial[0,keyPressFN:frameN+1] = 270 # down
        # registering response at the end of the trial for short trials:
        if key_arrow.status == STARTED and trialT <= 5 and t > trialT:
            theseKeys = event.getKeys(keyList=['left','right','up','down'])
            if len(theseKeys) > 0:
                dirFdbkL.setAutoDraw(True)
                dirFdbkR.setAutoDraw(True)
                someKeyPressed = True
                key_qn = True
                if ringDrawn:
                    respGiven = True
                if 'left' in theseKeys:
                    print '"left" key is pressed'
                    behRespTrial = 180
                    drawFdbkAngle(dirFdbkL, -1, behRespTrial)
                    drawFdbkAngle(dirFdbkR, 1, behRespTrial)
                elif 'right' in theseKeys:
                    print '"right" key is pressed'
                    behRespTrial = 0
                    drawFdbkAngle(dirFdbkL, -1, behRespTrial)
                    drawFdbkAngle(dirFdbkR, 1, behRespTrial)
                elif 'up' in theseKeys:
                    print '"up" key is pressed'
                    behRespTrial = 90
                    drawFdbkAngle(dirFdbkL, -1, behRespTrial)
                    drawFdbkAngle(dirFdbkR, 1, behRespTrial)
                elif 'down' in theseKeys:
                    print '"down" key is pressed'
                    behRespTrial = 270
                    drawFdbkAngle(dirFdbkL, -1, behRespTrial)
                    drawFdbkAngle(dirFdbkR, 1, behRespTrial)

        # after-trial question about the trial stability
        if not key_qn and t > trialT and stabQn:
            qntxtL.setAutoDraw(True)
            qntxtR.setAutoDraw(True)
            theseKeys = event.getKeys(keyList=['1','2','3','4'])
            if len(theseKeys)>0:
                if '1' in theseKeys:
                    qnResp = 1
                elif '2' in theseKeys:
                    qnResp = 2
                elif '3' in theseKeys:
                    qnResp = 3
                elif '4' in theseKeys:
                    qnResp = 4
                qntxtL.setAutoDraw(False)
                qntxtR.setAutoDraw(False)
                key_qn = True

        # after-trial question about the extent of the central motion pattern:
        if t > trialT and centTask:
            theseKeys = event.getKeys(keyList=['z','x','c'])
            if len(theseKeys)>0:
                ringL.setAutoDraw(True)
                ringR.setAutoDraw(True)
                szRelMax = np.max([szRelL,szRelR])
                if 'z' in theseKeys:
                    ringL = ringSzFn(ringL, periGap*szRelMax, 1) # increase the ring size
                    ringR = ringSzFn(ringR, periGap*szRelMax, 1)
                elif 'x' in theseKeys:
                    ringL = ringSzFn(ringL, periGap*szRelMax, -1) # decrease the ring size
                    ringR = ringSzFn(ringR, periGap*szRelMax, -1)
                elif 'c' in theseKeys:
                    ringL.setAutoDraw(False)
                    ringR.setAutoDraw(False)
                    ringSz = int(np.max([ringL.size[0],ringR.size[0]]))
                    print 'ringSz = ' + str(ringSz)
                    ringDrawn = True
                    if key_qn: # if the arrows are already drawn
                        respGiven = True

        if t > trialT and not elStopped:
            # stopping eye-tracking recording:
            if et:
                elEndRec(el)
                elStopped = True

        # pause text and data exporting
        if respGiven and not key_pause and t>trialT:
            if not behRespRecorded: # a flag for data recording
                # Make sure to record the release of a key at trial end
                if someKeyPressed and not centTask and trialT>5:
                    if whichKeyPressed == 'left':
                        behRespTrial[0,keyPressFN:(trialT*nFrames)] = 180
                    if whichKeyPressed == 'right':
                        behRespTrial[0,keyPressFN:(trialT*nFrames)] = 0
                    if whichKeyPressed == 'up':
                        behRespTrial[0,keyPressFN:(trialT*nFrames)] = 90
                    if whichKeyPressed == 'down':
                        behRespTrial[0,keyPressFN:(trialT*nFrames)] = 270
                    dirFdbkL.setAutoDraw(False)
                    dirFdbkR.setAutoDraw(False)
                    print '...released after ' + \
                        str(np.around(((trialT*nFrames)-keyPressFN)/60,2)) + 's'
                    print 'recorded post-trial response'
                # Recording the responses:
                pauseTextL.setAutoDraw(True)
                pauseTextR.setAutoDraw(True)
            if 'space' in event.getKeys(keyList=['space']):
                behRespRecorded = True
                # Computing and recording predominance:
                if trialT > 5:
                    nNa = np.count_nonzero(np.isnan(behRespTrial))
                    nf000 = np.count_nonzero(behRespTrial==0)
                    nf090 = np.count_nonzero(behRespTrial==90)
                    nf180 = np.count_nonzero(behRespTrial==180)
                    nf270 = np.count_nonzero(behRespTrial==270)
                    #print nNa, nf000, nf090, nf180, nf270
                else:
                    nNa = 0
                    nf000 = 0
                    nf090 = 0
                    nf180 = 0
                    nf270 = 0
                    if behRespTrial==0: nf000 = 1
                    if behRespTrial==90: nf090 = 1
                    if behRespTrial==180: nf180 = 1
                    if behRespTrial==270: nf270 = 1
                dT = pd.DataFrame({'expName': expName,
                                'time': expInfo['time'],
                                'participant': expInfo['participant'],
                                'session': expInfo['session'],
                                'trialN': nDone,
                                'dirL': dirL, 'dirR': dirR,
                                'vL': vL, 'vR': vR, 'szL': szL, 'szR': szR,
                                'sfL': sfL, 'sfR': sfR, 'BvL': BvL, 'BvR': BvR,
                                'BsfL': BsfL, 'BsfR': BsfR,
                                'colorL': str(colorL), 'colorR': str(colorR), 'sat': sat,
                                'fovGap': fovGap, 'fovFade': fovFade,
                                'periGap': periGap, 'periFade': periFade,
                                'szRelL': szRelL, 'szRelR': szRelR,
                                'offX': offX, 'offY': offY, 'tOffL': tOffL, 'tOffR': tOffR,
                                'trialT': trialT, 'nFrames': nFrames, 'nNa': nNa,
                                'nf000': nf000, 'nf090': nf090, 'nf180': nf180, 'nf270': nf270,
                                'pd000': [nf000 / (trialT * nFrames)],
                                'pd090': [nf090 / (trialT * nFrames)],
                                'pd180': [nf180 / (trialT * nFrames)],
                                'pd270': [nf270 / (trialT * nFrames)],
                                'qnResp': qnResp, 'ringSz': ringSz})
                # to preserve the column order:
                dataCols = ['expName', 'time', 'participant', 'session', 'trialN', 'dirL', 'dirR',
                            'vL', 'vR', 'szL', 'szR', 'sfL', 'sfR', 'tfL', 'tfR', 'BvL', 'BvR',
                            'BsfL', 'BsfR', 'colorL', 'colorR', 'sat', 'fovGap', 'fovFade', 
                            'periGap', 'periFade', 'szRelL', 'szRelR', 'offX', 'offY', 'tOffL', 'tOffR',
                            'trialT', 'nFrames', 'nNa', 'nf000', 'nf090', 'nf180', 'nf270', 
                            'pd000', 'pd090', 'pd180', 'pd270', 'qnResp', 'ringSz']
                if nDone == 1:
                    df = dT
                else:
                    df = pd.concat([df,dT])
                # Recording the data to a csv file:
                df.to_csv(dataFileName, index=False, columns=dataCols)
                print 'wrote the data set to ' + dataFileName
                print 'spacebar pressed - continuing to the next trial'
                pauseTextL.setAutoDraw(False)
                pauseTextR.setAutoDraw(False)
                key_pause = True

        # *ISI* period
        if ISI.status == NOT_STARTED and t>=trialT and key_pause:
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

# trialsFilePath = filePath + os.sep + fileName + '_trials'
# trials.saveAsPickle(trialsFilePath)
# trials.saveAsText(trialsFilePath)
# print trials

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
