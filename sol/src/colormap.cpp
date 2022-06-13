
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
#include <vtkNew.h>
#include <vtkColorSeries.h>
#include <vtkPolyDataMapper.h>
#include <vtkColorTransferFunction.h>
#include <vtkCommand.h>
#include <vtkAlgorithmOutput.h>
#include <vtkProgrammableFilter.h>



#pragma warning(pop)

const double red[]{ 247, 22, 22, 1 };
const double blue[]{ 20, 43, 222, 1 };


struct filter_params
{
	vtkDoubleArray* scalars;
	vtkPlaneSource* plane;
	vtkPolyDataMapper* mapper;
};

void timercallback(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))
{
	//vtkSmartPointer<vtkDoubleArray> scalars = static_cast<vtkDoubleArray*>(clientData);

	filter_params& p = *static_cast<filter_params*>(clientData);
	vtkSmartPointer<vtkRenderWindowInteractor> iren = static_cast<vtkRenderWindowInteractor*>(caller);

	iren->Render();
	p.scalars->Modified();
	//p.mapper->Update();
	//iren->Render();
}



namespace LookupTables
{
	auto RedBlue()
	{
		// Map the scalar values in the image to colors with a lookup table:
		vtkNew<vtkLookupTable> lookupTable;
		//lookupTable->SetRange(-0.5, 0.5);
		lookupTable->SetHueRange(.05, .7);
		lookupTable->SetSaturationRange(.7, .7);
		lookupTable->SetValueRange(.6, .6);
		lookupTable->SetRampToLinear();

		//lookupTable->SetNumberOfTableValues(1024);
		lookupTable->SetTableRange(-1., 1.);
		//lookupTable->SetTableValue(0, red);
		//lookupTable->SetTableValue(lookupTable->GetNumberOfTableValues() - 1, blue);// lookupTable->GetNumberOfTableValues() -
		lookupTable->IndexedLookupOff();
		lookupTable->Build();

		return lookupTable;
	}
};

void update_filter(void* arguments)
{
	filter_params &p = *static_cast<filter_params*>(arguments);
	//p.filter->GetPolyDataOutput()->GetPointData()->SetScalars(p.data);
	//p.filter->GetPolyDataOutput()->(p.filter->GetInput());
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
		plane->SetXResolution(dims[0] - 1);
		plane->SetYResolution(dims[1] - 1);
		plane->Update();
		plane->GetOutput()->GetPointData()->SetScalars(scalars);


		vtkNew<vtkLookupTable> lookupTable;
		vtkNew<vtkColorSeries> colorSeries;
		colorSeries->SetColorScheme(vtkColorSeries::ColorSchemes::CITRUS);
		colorSeries->BuildLookupTable(lookupTable, vtkColorSeries::ORDINAL);
		

		//vtkNew<vtkSimpleElevationFilter> filter;
		//filter->SetInputConnection(image->GetOutputPort());
		//filter->SetVector(1, 1, 1);

		//scalarValuesToColors->SetLookupTable(lookupTable);
		//scalarValuesToColors->PassAlphaToOutputOn();
		//scalarValuesToColors->SetInputData(image);

		vtkNew<vtkPolyDataMapper> mapper;
		mapper->SetInputConnection(plane->GetOutputPort());
		mapper->SetScalarRange(-.5, .5);
		//mapper->ScalarVisibilityOn();
		mapper->SetLookupTable(lookupTable);
		//plotMapper->GetProperty()->SetInterpolationTypeToNearest();

		vtkNew<vtkActor> plotActor;
		plotActor->SetMapper(mapper);

		vtkNew<vtkRenderer> renderer;
		vtkNew<vtkRenderWindow> renderWindow;
		vtkNew<vtkInteractorStyleImage> style;
		vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;

		renderWindow->AddRenderer(renderer);
		renderWindowInteractor->SetInteractorStyle(style);
		renderWindowInteractor->SetRenderWindow(renderWindow);
		
		renderWindow->SetSize(500, 500);
		renderWindow->SetWindowName("SymPhas");


		// first initialize the render interactor
		renderWindowInteractor->Initialize();

		filter_params* p = new filter_params{ scalars, plane, mapper };

		// now we can add a timer based on the callback method
		vtkNew<vtkCallbackCommand> timerCallback;
		timerCallback->SetCallback(timercallback);
		timerCallback->SetClientData(p);

		renderWindowInteractor->AddObserver(vtkCommand::TimerEvent, timerCallback);
		renderWindowInteractor->CreateRepeatingTimer(std::max(params::viz_interval, 20));

		// Visualize
		renderer->AddActor(plotActor);
		renderer->ResetCamera();
		

		renderWindow->Render();
		renderWindowInteractor->Start();
	}
}



