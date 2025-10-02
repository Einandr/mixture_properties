import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import os
from enum import Enum
from dataclasses import dataclass

from . import component_terra as tc
from . import component_chemkin as cc
from . import component_mixture as mc
from . import output as ou

# Константы:
# R0 = 8.31

class Source(Enum):
    T = 0    # TERRA
    C = 1    # CHEMKIN

@dataclass
class Component:
    value: float
    source: Source


class Material:
    def __init__(self, mixture_reference, output_file_name, df_terra, components_chemkin, show_plots, is_gas, T0, T_base, dT_phase_transition, T_first, T_last, dT):
        self.mixture_reference = mixture_reference
        self.output_file_name = output_file_name
        self.T_grid_mixture_nonunique = []
        self.components = components_chemkin        # Объединенный список компонентов CHEMKIN и TERRA
        self.component_terra_names_list = []
        self.df_terra = df_terra
        # self.components_chemkin = components_chemkin
        self.show_plots = show_plots
        self.is_gas = is_gas
        self.T0 = T0
        self.T_base = T_base
        self.dT_phase_transition = dT_phase_transition
        self.T_first = T_first
        self.T_last = T_last
        self.dT = dT
        self.path_terra = 'terra'
        os.makedirs(self.path_terra, exist_ok=True)
        os.chdir(self.path_terra)
        for key, value in self.mixture_reference.items():
            for key2, in_component in value.items():
                print('inside material key2 is:', key2)
                print('inside material in_component is:', in_component)
                if in_component.source == Source.T:
                    component = tc.initialize_terra_component(key2, self.df_terra, self.T0, self.T_base, self.dT_phase_transition, self.T_first, self.T_last, self.dT)
                    if component.name not in self.component_terra_names_list:
                        print('\nadding component with name ', component.name, '\n')
                        self.components.append(component)
                        self.component_terra_names_list.append(component.name)
                        self.T_grid_mixture_nonunique.extend(tc.show_properties_of_component(component, self.show_plots,  self.T_base))
                # elif in_component.source == Source.C:
                #     component = cc.initialize_chemkin_component(key2, self.df_terra, self.T0, self.T_base, self.dT_phase_transition, self.T_first, self.T_last, self.dT)




        print('CLASS material components names list:', self.component_terra_names_list)

        for c in self.components:
            print(c.name)
        self.T_grid = pd.Series(self.T_grid_mixture_nonunique).drop_duplicates().sort_values(ascending=True).tolist()
        # for objects of class Terra Component
        self.vector_get_mixture_Cp = np.vectorize(self.get_mixture_Cp, excluded=['mass_fractions_dict'])
        self.vector_get_mixture_H = np.vectorize(self.get_mixture_H, excluded=['mass_fractions_dict'])
        self.vector_get_mixture_viscosity = np.vectorize(self.get_mixture_viscosity, excluded=['mass_fractions_dict'])
        self.vector_get_mixture_heat_conductivity = np.vectorize(self.get_mixture_heat_conductivity, excluded=['mass_fractions_dict'])
        self.vector_get_mixture_diffusivity = np.vectorize(self.get_mixture_diffusivity, excluded=['mass_fractions_dict'])

        # for objects of class Mixture Component
        self.vector_get_heat_capacity = np.vectorize(self.get_heat_capacity, excluded=['mass_fractions'])
        self.vector_get_enthalpy = np.vectorize(self.get_enthalpy, excluded=['mass_fractions'])
        self.vector_get_viscosity = np.vectorize(self.get_viscosity, excluded=['mass_fractions'])
        self.vector_get_heat_conductivity = np.vectorize(self.get_heat_conductivity, excluded=['mass_fractions'])
        self.vector_get_diffusivity = np.vectorize(self.get_diffusivity, excluded=['mass_fractions'])

        os.chdir('..')
        self.path_mixture = 'mixture'
        os.makedirs(self.path_mixture, exist_ok=True)
        os.chdir(self.path_mixture)

        self.writer = pd.ExcelWriter(''.join((self.output_file_name, '.xlsx')), engine="xlsxwriter")
        self.list_of_constant_properties = []
        self.dict_of_constant_properties = {}
        self.dict_of_T_dependent_properties = {}
        for mixture_name, mass_fractions_dict in self.mixture_reference.items():
            # data_constant_props = calculate_mixture_constant_properties(name, mixture)
            data_constant_props = self.__calculate_mixture_properties_constant(mixture_name, mass_fractions_dict)
            self.list_of_constant_properties.append(data_constant_props)
            self.dict_of_constant_properties[mixture_name] = data_constant_props
            # data_T_dependent_properties = calculate_mixture_properties(writer_gas, name, mixture)
            data_T_dependent_properties = self.__calculate_mixture_properties_Tdependent(mixture_name, mass_fractions_dict)
            self.dict_of_T_dependent_properties[mixture_name] = data_T_dependent_properties
        self.constant_props = pd.concat(self.list_of_constant_properties)
        self.constant_props.to_excel(self.writer, index=True, sheet_name='constant_props')
        self.writer.close()



        self.dict_of_components = {}
        for mixture_name, mass_fractions_dict in self.mixture_reference.items():
            mu = self.dict_of_constant_properties[mixture_name]['mu'].iloc[0]
            rho = self.dict_of_constant_properties[mixture_name]['rho'].iloc[0]
            dH0 = self.dict_of_constant_properties[mixture_name]['dH0'].iloc[0]
            T_grid = self.dict_of_T_dependent_properties[mixture_name]['T [K]'].to_numpy()
            Cp_grid = self.dict_of_T_dependent_properties[mixture_name]['Cp [Дж/кг-К]'].to_numpy()
            H_grid = self.dict_of_T_dependent_properties[mixture_name]['H [Дж/кг]'].to_numpy()
            if not self.is_gas:
                mixture_component = mc.initialize_mixture_component_dispersed(mixture_reference, mixture_name, self.T_base, mu, rho, dH0, T_grid, Cp_grid, H_grid)
                self.dict_of_components[mixture_name] = mixture_component
            else:
                viscosity_grid = self.dict_of_T_dependent_properties[mixture_name]['Mu_visc [kg/m-s]'].to_numpy()
                heat_conductivity_grid = self.dict_of_T_dependent_properties[mixture_name]['Lambda [W/m-K]'].to_numpy()
                diffusivity_grid = self.dict_of_T_dependent_properties[mixture_name]['D [m^2/s]'].to_numpy()
                mixture_component = mc.initialize_mixture_component_gas(mixture_reference, mixture_name, self.T_base, mu, rho, dH0, T_grid, Cp_grid, H_grid, viscosity_grid, heat_conductivity_grid, diffusivity_grid)
                self.dict_of_components[mixture_name] = mixture_component


        os.chdir('..')

    # МЕТОДЫ будут вызывать в зависимости от типа компонента - терра либо чемкин, соответствующие внутренние методы компонента с одинаковыми названиями
    def get_mixture_Cp(self, mass_fractions_dict, T):
        result = 0
        for key, component in mass_fractions_dict.items():
            for comp in self.components:
                try:
                    if key == comp.name:
                        result += comp.get_Cp(T) * component.value
                except ValueError:
                    print('no name for component ', comp)
                continue
        return result


    def get_mixture_viscosity(self, mass_fractions_dict, T):
        result = 0
        for key, component in mass_fractions_dict.items():
            for comp in self.components:
                if key == comp.name:
                    result += comp.get_viscosity(T) * component.value
        return result

    def get_mixture_heat_conductivity(self, mass_fractions_dict, T):
        result = 0
        for key, component in mass_fractions_dict.items():
            for comp in self.components:
                if key == comp.name:
                    result += comp.get_heat_conductivity(T) * component.value
        return result

    def get_mixture_diffusivity(self, mass_fractions_dict, T):
        result = 0
        for key, component in mass_fractions_dict.items():
            for comp in self.components:
                if key == comp.name:
                    result += comp.get_diffusivity(T) * component.value
        return result


    def get_mixture_H(self, mass_fractions_dict, T):
        result = 0
        for key, component in mass_fractions_dict.items():
            for comp in self.components:
                if key == comp.name:
                    result += comp.get_H(T) * component.value
        return result

    def get_mixture_mu(self, mass_fractions_dict):
        result = 0
        for key, component in mass_fractions_dict.items():
            for comp in self.components:
                if key == comp.name:
                    # print('mixture components mu: ', key, comp.get_mu(), value)
                    result += component.value / comp.get_mu()
        return 1 / result

    def get_mixture_rho(self, mass_fractions_dict):
        result = 0
        # prevent division by zero
        if not self.is_gas:
            for key, component in mass_fractions_dict.items():
                for comp in self.components:
                    if key == comp.name:
                        result += component.value / comp.get_rho()
            return 1 / result
        else:
            return 0

    def get_mixture_dH0(self, mass_fractions_dict):
        result = 0
        for key, component in mass_fractions_dict.items():
            for comp in self.components:
                if key == comp.name:
                    result += comp.get_dH0() * component.value
        return result

    def __calculate_mixture_properties_Tdependent(self, mixture_name, mass_fractions_dict):
        Cp_grid = self.vector_get_mixture_Cp(mass_fractions_dict, self.T_grid)
        H_grid = self.vector_get_mixture_H(mass_fractions_dict, self.T_grid)
        if self.is_gas:
            viscosity_grid = self.vector_get_mixture_viscosity(mass_fractions_dict, self.T_grid)
            heat_conductivity_grid = self.vector_get_mixture_heat_conductivity(mass_fractions_dict, self.T_grid)
            diffusivity_grid = self.vector_get_mixture_diffusivity(mass_fractions_dict, self.T_grid)

        if self.show_plots:
            # plots
            if not self.is_gas:
                plt.figure(1)
                plt.grid(True)
                plt.plot(self.T_grid, Cp_grid, label=r'$Ср$')
                plt.legend(fontsize='medium', loc='best')
                plt.title(r'$Теплоёмкость$')
                plt.xlabel(r'$Т\ (К)$')
                plt.ylabel(r'$C_p\ (\frac{Дж}{кг\cdot К})$')
                plt.xlim(0)
                plt.ylim(0)
                plt.savefig(''.join(('disp_', mixture_name, '_Cp.jpeg')), dpi=400, bbox_inches='tight')
                plt.close()

                plt.figure(2)
                plt.grid(True)
                plt.plot(self.T_grid, H_grid, label=r'$H$')
                plt.legend(fontsize='medium', loc='best')
                plt.title(r'$Энтальпия$')
                plt.xlabel(r'$Т\ (К)$')
                plt.ylabel(r'$H\ (\frac{Дж}{кг})$')
                plt.xlim(0)
                plt.ylim(0)
                plt.savefig(''.join(('disp_', mixture_name, '_H.jpeg')), dpi=400, bbox_inches='tight')
                plt.close()
            else:
                plt.figure(1)
                plt.grid(True)
                plt.plot(self.T_grid, Cp_grid, label=r'$Ср$')
                plt.legend(fontsize='medium', loc='best')
                plt.title(r'$Теплоёмкость$')
                plt.xlabel(r'$Т\ (К)$')
                plt.ylabel(r'$C_p\ (\frac{Дж}{кг\cdot К})$')
                plt.xlim(0)
                plt.ylim(0)
                plt.savefig(''.join(('gas_', mixture_name, '_Cp.jpeg')), dpi=400, bbox_inches='tight')
                plt.close()

                plt.figure(2)
                plt.grid(True)
                plt.plot(self.T_grid, H_grid, label=r'$H$')
                plt.legend(fontsize='medium', loc='best')
                plt.title(r'$Энтальпия$')
                plt.xlabel(r'$Т\ (К)$')
                plt.ylabel(r'$H\ (\frac{Дж}{кг})$')
                plt.xlim(0)
                plt.ylim(0)
                plt.savefig(''.join(('gas_', mixture_name, '_H.jpeg')), dpi=400, bbox_inches='tight')
                plt.close()

                plt.figure(3)
                plt.grid(True)
                plt.plot(self.T_grid, viscosity_grid, label=r'$\mu$')
                plt.legend(fontsize='medium', loc='best')
                plt.title(r'$Вязкость$')
                plt.xlabel(r'$Т\ (К)$')
                plt.ylabel(r'$\mu\ (Па \cdot с)$')
                plt.xlim(0)
                plt.ylim(0)
                plt.savefig(''.join(('gas_', mixture_name, '_viscosity.jpeg')), dpi=400, bbox_inches='tight')
                plt.close()

                plt.figure(4)
                plt.grid(True)
                plt.plot(self.T_grid, heat_conductivity_grid, label=r'$\lambda$')
                plt.legend(fontsize='medium', loc='best')
                plt.title(r'$Теплопроводность$')
                plt.xlabel(r'$Т\ (К)$')
                plt.ylabel(r'$\lambda\ (\frac{Вт}{м \cdot K})$')
                plt.xlim(0)
                plt.ylim(0)
                plt.savefig(''.join(('gas_', mixture_name, '_heat_conductivity.jpeg')), dpi=400, bbox_inches='tight')
                plt.close()

                plt.figure(5)
                plt.grid(True)
                plt.plot(self.T_grid, diffusivity_grid, label=r'$D$')
                plt.legend(fontsize='medium', loc='best')
                plt.title(r'$Диффузия$')
                plt.xlabel(r'$Т\ (К)$')
                plt.ylabel(r'$D\ (\frac{м^2}{c})$')
                plt.xlim(0)
                plt.ylim(0)
                plt.savefig(''.join(('gas_', mixture_name, '_diffusivity.jpeg')), dpi=400, bbox_inches='tight')
                plt.close()


        if not self.is_gas:
            frame = {'T [K]': pd.Series(self.T_grid),
                     'Cp [Дж/кг-К]': pd.Series(Cp_grid),
                     'H [Дж/кг]': pd.Series(H_grid)}
        else:
            frame = {'T [K]': pd.Series(self.T_grid),
                     'Cp [Дж/кг-К]': pd.Series(Cp_grid),
                     'H [Дж/кг]': pd.Series(H_grid),
                     'Mu_visc [kg/m-s]': pd.Series(viscosity_grid),
                     'Lambda [W/m-K]': pd.Series(heat_conductivity_grid),
                     'D [m^2/s]': pd.Series(diffusivity_grid)}


        # data = pd.DataFrame(frame)
        data = pd.DataFrame(frame)
        data.set_index('T [K]', drop=False, inplace=True)
        print('data check: ', data)
        data.to_excel(self.writer, index=False, sheet_name=mixture_name)
        ou.output_props(mixture_name, self.is_gas, self.T_base, data)
        return data

    def __calculate_mixture_properties_constant(self, mixture_name, mass_fractions_dict):
        mu = self.get_mixture_mu(mass_fractions_dict)
        rho = self.get_mixture_rho(mass_fractions_dict)
        dH0 = self.get_mixture_dH0(mass_fractions_dict)
        return pd.DataFrame(data=[[mu, rho, dH0]], index=[mixture_name], columns=['mu', 'rho', 'dH0'])

    # next methods no need for check in self.dict_of_components - this is for MIXTURE COMPONENTS
    def get_density(self, mass_fractions):
        result = 0
        for mf_name, value in mass_fractions.items():
            result += value / self.dict_of_components[mf_name].get_rho()
        return 1 / result

    def get_heat_capacity(self, mass_fractions, T):
        result = 0
        for mf_name, value in mass_fractions.items():
            result += value * self.dict_of_components[mf_name].get_Cp(T)
        return result

    def get_enthalpy(self, mass_fractions, T):
        result = 0
        for mf_name, value in mass_fractions.items():
            result += value * self.dict_of_components[mf_name].get_H(T)
        # print('get enthalpy result: ', result)
        return result

    def get_viscosity(self, mass_fractions, T):
        result = 0
        for mf_name, value in mass_fractions.items():
            result += value * self.dict_of_components[mf_name].get_viscosity(T)
        return result

    def get_heat_conductivity(self, mass_fractions, T):
        result = 0
        for mf_name, value in mass_fractions.items():
            result += value * self.dict_of_components[mf_name].get_heat_conductivity(T)
        return result

    def get_diffusivity(self, mass_fractions, T):
        result = 0
        for mf_name, value in mass_fractions.items():
            result += value * self.dict_of_components[mf_name].get_diffusivity(T)
        return result

    def get_formation_heat(self, mass_fractions):
        result = 0
        for mf_name, value in mass_fractions.items():
            result += value * self.dict_of_components[mf_name].get_dH0()
        return result

    def get_mu(self, mass_fractions):
        denom = 0
        for mf_name, value in mass_fractions.items():
            mu_mixture_component = self.dict_of_components[mf_name].get_mu()
            # print('mu of mixture component: mf_name, value, mu: ', mf_name, value, mu_mixture_component)
            denom += value/mu_mixture_component
        result = 1/denom
        return result

    def interpolated_value(self, x, x0, x1, y0, y1):
        return y0 + (y1 - y0) / (x1 - x0) * (x - x0)

    def get_temperature_from_enthalpy(self, mass_fractions, enthalpy):
        Cp_grid_mixture = self.vector_get_heat_capacity(mass_fractions, self.T_grid)
        H_grid_mixture = self.vector_get_enthalpy(mass_fractions, self.T_grid)

        # print('T_grid: ', self.T_grid)
        # print('Cp_grid_mixture: ', Cp_grid_mixture)
        # print('H_grid_mixture: ', H_grid_mixture)
        max_index = len(self.T_grid) - 1
        H = enthalpy
        if H <= H_grid_mixture[0]:
            # T = H/Cp
            temperature = H / Cp_grid_mixture[0]
            # print('reconstructing temperature from enthalpy - under lower table range')
        elif H >= H_grid_mixture[max_index]:
            # dH = H-H_back
            # dH = dT * Cp_back
            # T-T_back = H-H_back / Cp_back
            # T = T_back + (H-H_back / Cp_back)
            temperature = self.T_grid[max_index] + (H-H_grid_mixture[max_index]) / Cp_grid_mixture[max_index]
            # print('reconstructing temperature from enthalpy - over upper table range', temperature, enthalpy)
        else:
            it = np.searchsorted(H_grid_mixture, H, side='left')
            T_0 = self.T_grid[it - 1]
            T_1 = self.T_grid[it]
            H_0 = H_grid_mixture[it - 1]
            H_1 = H_grid_mixture[it]
            temperature = self.interpolated_value(H, H_0, H_1, T_0, T_1)
            # print('reconstructing temperature from enthalpy, iteration, T_0, T_1, H_0, H_1: ', it, T_0, T_1, H_0, H_1)
            print('components in interpolation t from enthalpy: T1 T2 temperature H1 H2 H: ', T_0, T_1, temperature, H_0, H_1, H)

        return temperature


