
/* ***************************************************************************
 * This file is part of the SymPhas library, a framework for implementing
 * solvers for phase-field problems with compile-time symbolic algebra.
 * 
 * Copyright (c) 2018-2021 by Steven A. Silber and Mikko Karttunen
 * 
 * SymPhas is free software, which can be redistributed or modified under
 * the terms of the GNU Lesser General Public License (LGPL) as published
 * by the Free Software Foundation; LGPL version 3, or later versions at
 * your choice.
 * 
 * SymPhas is distributed with the faith that it will be helpful and
 * practical but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
 * General Public License for more details.
 *
 * ***************************************************************************
 */


#include "colormap.h"
#pragma once
#pragma warning(push, 0)   


#include <vtkInteractorStyleImage.h>
#include <vtkLookupTable.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkUniformGrid.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkProgrammableFilter.h>
#include <vtkCallbackCommand.h>
#include <vtkDataSetMapper.h>
#include <vtkSimpleElevationFilter.h>
#include <vtkPlaneSource.h>

#pragma warning(pop)


void timercallback(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))
{
	vtkSmartPointer<vtkInteractorStyleImage> programmableFilter = static_cast<vtkInteractorStyleImage*>(clientData);
	vtkRenderWindowInteractor* iren = static_cast<vtkRenderWindowInteractor*>(caller);

	programmableFilter->Modified();
	iren->Render();
}

void ColourPlot2d::init(scalar_t* values, len_type* dims)
{
	
	params::viz_interval = 4;
	if (params::viz_interval > 0)
	{
		size_t len = grid::length<2>(dims);

		vtkNew<vtkDoubleArray> scalars;
		scalars->SetArray(values, len, 1);

		vtkNew<vtkPlaneSource> plane;
		plane->SetXResolution(dims[0]);
		plane->SetYResolution(dims[1]);
		plane->Update();
		plane->GetOutput()->GetPointData()->SetScalars(scalars);


		// Map the scalar values in the image to colors with a lookup table:
		vtkNew<vtkLookupTable> lookupTable;
		lookupTable->SetNumberOfTableValues(1024);
		lookupTable->SetRange(-0.5, 0.5);
		lookupTable->Build();

		//vtkNew<vtkSimpleElevationFilter> filter;
		//filter->SetInputConnection(image->GetOutputPort());
		//filter->SetVector(1, 1, 1);

		//scalarValuesToColors->SetLookupTable(lookupTable);
		//scalarValuesToColors->PassAlphaToOutputOn();
		//scalarValuesToColors->SetInputData(image);

		// Create an image actor
		vtkNew<vtkDataSetMapper> magnitude;
		magnitude->SetInputConnection(plane->GetOutputPort());
		magnitude->SetScalarRange(-1, 1);
		magnitude->ScalarVisibilityOn();
		magnitude->SetLookupTable(lookupTable);
		//plotMapper->GetProperty()->SetInterpolationTypeToNearest();

		vtkNew<vtkActor> plotActor;
		plotActor->SetMapper(magnitude);

		vtkNew<vtkRenderer> renderer;
		vtkNew<vtkRenderWindow> renderWindow;
		vtkNew<vtkInteractorStyleImage> style;
		vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;

		renderWindow->AddRenderer(renderer);
		//renderWindowInteractor->SetInteractorStyle(style);
		renderWindowInteractor->SetRenderWindow(renderWindow);

		// first initialize the render interactor
		///renderWindowInteractor->Initialize();

		// now we can add a timer based on the callback method
		//renderWindowInteractor->CreateRepeatingTimer(std::max(params::viz_interval, 50));
		//vtkNew<vtkCallbackCommand> timerCallback;
		//timerCallback->SetCallback(timercallback);
		//timerCallback->SetClientData(scalarValuesToColors);

		//renderWindowInteractor->AddObserver(vtkCommand::TimerEvent, timerCallback);

		// Visualize
		renderer->AddActor(plotActor);
		renderer->ResetCamera();
		
		

		renderWindow->Render();
		renderWindowInteractor->Start();
	}
}



