#pvpython
from paraview.simple import *
import csv

paraview.simple._DisableFirstRenderCameraReset()

solutionpvd = PVDReader(FileName="results/solution.pvd")
SetActiveSource(solutionpvd)

tsteps = solutionpvd.TimestepValues
line = PlotOverLine( Source="High Resolution Line Source")
DataRepresentation7 = Show()

line.Source.Point1 = [0, 5e-7, 0.0]
line.Source.Point2 = [0.0005000000237487257, 5e-7, 0.0]
line.Source.Resolution = 5000 #Uniform sampling
component = 0
sample_size = 0.5e-3
delta_h = sample_size/line.Source.Resolution

with open('length5000.csv', 'w') as f:
    writer = csv.writer(f)
    for TimeStepNum in range(0,len(tsteps)):
        view = GetActiveView()
        view.ViewTime = tsteps[TimeStepNum]
        Render()
        fetchData = paraview.servermanager.Fetch(line)
        pointData = fetchData.GetPointData()
        fieldData = pointData.GetArray("Damage: $s$")
        dam = [fieldData.GetComponent(i,component) for i in range(fieldData.GetSize())]
        filtered_dam = [x for x in dam if x < 0.5]
        writer.writerow([tsteps[TimeStepNum], (len(filtered_dam) * delta_h)])
