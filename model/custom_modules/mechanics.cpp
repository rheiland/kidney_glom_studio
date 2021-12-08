#include "./mechanics.h"

void epithelial_special_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	// BM adhesion 
		// is it time to detach (attachment lifetime)
		// am I unattached by capable? 
			// search through neighbors, find closest BM type agent 
			// form adhesion 
		// elastic adhesion 
	
	// plasto-elastic. 
		// elastic: movement towards rest position 
	
	static int nRP = 0; // "rest_position"
	
	std::vector<double> displacement = pCell->custom_data.vector_variables[nRP].value ; 
	
	
	
	
	
	
	
	
	
	return; 
}


void plasto_elastic_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	// BM adhesion 
		// is it time to detach (attachment lifetime)
		// am I unattached by capable? 
			// search through neighbors, find closest BM type agent 
			// form adhesion 
		// elastic adhesion 
	
	// plasto-elastic. 
		// elastic: movement towards rest position 
	
	static int nRP = 0; // "rest_position"
	
	// displacement 
	std::vector<double> displacement = pCell->custom_data.vector_variables[nRP].value - pCell->position; 
	
	static int nEConst = pCell->custom_data.find_variable_index( "cell_elasticity" );
	static int nPConst = pCell->custom_data.find_variable_index( "cell_plasticity" );

	// first, update the agent's velocity based upon the elastic model
	axpy( &( pCell->velocity ) , pCell->custom_data[nEConst] , displacement );

	// now, plastic mechanical relaxation

	double plastic_temp_constant = -dt * pCell->custom_data[nPConst];
	axpy( &(pCell->custom_data.vector_variables[nRP].value) , plastic_temp_constant , displacement );
	
	return; 
}

// specialized potential function 

void add_heterotypic_potentials(Cell* my_cell , Cell* other_agent)
{
	static int nCellID = my_cell->custom_data.find_variable_index( "cell_ID" ); 
	
	// if( this->ID == other_agent->ID )
	if( my_cell == other_agent )
	{ return; }
	
	static int nRAOC = my_cell->custom_data.find_variable_index( "relative_adhesion_other_cells" ); 
	static int nRAOCT = my_cell->custom_data.find_variable_index( "relative_adhesion_other_cell_types" ); 
	
	double rel_heterotypic_adhesion = my_cell->custom_data[nRAOCT];  
	double rel_other_cells_adhesion = my_cell->custom_data[nRAOC];  
	
	// 12 uniform neighbors at a close packing distance, after dividing out all constants
	static double simple_pressure_scale = 0.027288820670331; // 12 * (1 - sqrt(pi/(2*sqrt(3))))^2 
	// 9.820170012151277; // 12 * ( 1 - sqrt(2*pi/sqrt(3)))^2

	double distance = 0; 
	for( int i = 0 ; i < 3 ; i++ ) 
	{ 
		my_cell->displacement[i] = my_cell->position[i] - (*other_agent).position[i]; 
		distance += (my_cell->displacement[i]) * (my_cell->displacement[i]); 
	}
	// Make sure that the distance is not zero
	
	distance = std::max(sqrt(distance), 0.00001); 
	
	//Repulsive
	double R = my_cell->phenotype.geometry.radius + (*other_agent).phenotype.geometry.radius; 
	
	double RN = my_cell->phenotype.geometry.nuclear_radius + (*other_agent).phenotype.geometry.nuclear_radius;	
	double temp_r, c;
	if( distance > R ) 
	{
		temp_r=0;
	}
	else
	{
		temp_r = -distance; // -d
		temp_r /= R; // -d/R
		temp_r += 1.0; // 1-d/R
		temp_r *= temp_r; // (1-d/R)^2 
		
		// add the relative pressure contribution 
		my_cell->state.simple_pressure += ( temp_r / simple_pressure_scale ); // New July 2017 
	}
	
	// August 2017 - back to the original if both have same coefficient 
	double effective_repulsion = sqrt( my_cell->phenotype.mechanics.cell_cell_repulsion_strength * other_agent->phenotype.mechanics.cell_cell_repulsion_strength );
	temp_r *= effective_repulsion; 
	
	//////////////////////////////////////////////////////////////////
	
	// Adhesive
	double max_interactive_distance = my_cell->phenotype.mechanics.relative_maximum_adhesion_distance * my_cell->phenotype.geometry.radius + 
		(*other_agent).phenotype.mechanics.relative_maximum_adhesion_distance * (*other_agent).phenotype.geometry.radius;
		
	if(distance < max_interactive_distance ) 
	{	
		// double temp_a = 1 - distance/max_interactive_distance; 
		double temp_a = -distance; // -d
		temp_a /= max_interactive_distance; // -d/S
		temp_a += 1.0; // 1 - d/S 
		temp_a *= temp_a; // (1-d/S)^2 
		
		// August 2017 - back to the original if both have same coefficient 
		double effective_adhesion = sqrt( my_cell->phenotype.mechanics.cell_cell_adhesion_strength * other_agent->phenotype.mechanics.cell_cell_adhesion_strength ); 
		
		int my_id = (int) my_cell->custom_data[nCellID] ; 
		int other_id = (int) other_agent->custom_data[nCellID] ; 
	
		if( my_id != other_id )
		{ effective_adhesion *= rel_other_cells_adhesion; }
		
		if( my_cell->type != other_agent->type )
		{ effective_adhesion *= rel_heterotypic_adhesion; }
		temp_a *= effective_adhesion; 
		
		temp_r -= temp_a;
	}
	/////////////////////////////////////////////////////////////////
	if( fabs(temp_r) < 1e-16 )
	{ return; }
	temp_r /= distance;
	axpy( &(my_cell->velocity) , temp_r , my_cell->displacement ); 
	
	return;
}

void heterotypic_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt)
{
	if( pCell->functions.add_cell_basement_membrane_interactions )
	{
		pCell->functions.add_cell_basement_membrane_interactions(pCell, phenotype,dt);
	}
	
	pCell->state.simple_pressure = 0.0; 
	
	//First check the neighbors in my current voxel
	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for(neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{
		// pCell->add_potentials(*neighbor);
		add_heterotypic_potentials( pCell, *neighbor ); 
	}
	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();

	for( neighbor_voxel_index = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
		neighbor_voxel_index != neighbor_voxel_index_end; 
		++neighbor_voxel_index )
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
			continue;
		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
		{
			// pCell->add_potentials(*neighbor);
			add_heterotypic_potentials( pCell, *neighbor ); 
		}
	}

	pCell->update_motility_vector(dt); 
	pCell->velocity += phenotype.motility.motility_vector; 
	
	return; 
}



void BM_special_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	return; 
}

