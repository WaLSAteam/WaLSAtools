# -----------------------------------------------------------------------------------------------------
# WaLSAtools - Wave analysis tools
# Copyright (C) 2025 WaLSA Team - Shahin Jafarzadeh et al.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 
# Note: If you use WaLSAtools for research, please consider citing:
# Jafarzadeh, S., Jess, D. B., Stangalini, M. et al. 2025, Nature Reviews Methods Primers, 5, 21
# -----------------------------------------------------------------------------------------------------

import ipywidgets as widgets # type: ignore
from IPython.display import display, clear_output, HTML # type: ignore
import os # type: ignore

# Function to detect if it's running in a notebook or terminal
def is_notebook():
    try:
        # Check if IPython is in the environment and if we're using a Jupyter notebook
        from IPython import get_ipython # type: ignore
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True  # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # IPython terminal
        else:
            return False  # Other types
    except (NameError, ImportError):
        # NameError: get_ipython is not defined
        # ImportError: IPython is not installed
        return False  # Standard Python interpreter or shell

# Function to print logo and credits for both environments
def print_logo_and_credits():
    logo_terminal = r"""
        __          __          _          _____            
        \ \        / /         | |        / ____|     /\    
         \ \  /\  / /  ▄▄▄▄▄   | |       | (___      /  \   
          \ \/  \/ /   ▀▀▀▀██  | |        \___ \    / /\ \  
           \  /\  /   ▄██▀▀██  | |____    ____) |  / ____ \ 
            \/  \/    ▀██▄▄██  |______|  |_____/  /_/    \_\
    """
    credits_terminal = """
        © WaLSA Team (www.WaLSA.team)
        -----------------------------------------------------------------------
        WaLSAtools v1.0.1 - Wave analysis tools
        Documentation: www.WaLSA.tools
        GitHub repository: www.github.com/WaLSAteam/WaLSAtools
        -----------------------------------------------------------------------
        If you use WaLSAtools in your research, please cite:   
        Jafarzadeh, S., Jess, D. B., Stangalini, M. et al. 2025, Nature Reviews Methods Primers, 5, 21 

        Free access to a view-only version: https://WaLSA.tools/nrmp     
        Supplementary Information: https://WaLSA.tools/nrmp-si     
        -----------------------------------------------------------------------
        Choose a category, data type, and analysis method from the list below,
        to get hints on the calling sequence and parameters:
        """

    credits_notebook = """
        <div style="margin-left: 30px; margin-top: 20px; font-size: 1.1em; line-height: 0.8;">
            <p>© WaLSA Team (<a href="https://www.WaLSA.team" target="_blank">www.WaLSA.team</a>)</p>
            <hr style="width: 70%; margin: 0; border: 0.98px solid #888; margin-bottom: 10px;">
            <p><strong>WaLSAtools</strong> v1.0.1 - Wave analysis tools</p>
            <p>Documentation: <a href="https://www.WaLSA.tools" target="_blank">www.WaLSA.tools</a></p>
            <p>GitHub repository: <a href="https://www.github.com/WaLSAteam/WaLSAtools" target="_blank">www.github.com/WaLSAteam/WaLSAtools</a></p>
            <hr style="width: 70%; margin: 0; border: 0.98px solid #888; margin-bottom: 10px;">
            <p>If you use <strong>WaLSAtools</strong> in your research, please cite:</p>
            <p>Jafarzadeh, S., Jess, D. B., Stangalini, M. et al. 2025, <em>Nature Reviews Methods Primers</em>, 5, 21</p>
            <p>Free access to a view-only version: <a href="https://WaLSA.tools/nrmp" target="_blank">www.WaLSA.tools</a></p>
            <p>Supplementary Information: <a href="https://WaLSA.tools/nrmp-si" target="_blank">www.WaLSA.tools</a></p>
            <hr style="width: 70%; margin: 0; border: 0.98px solid #888; margin-bottom: 15px;">
            <p>Choose a category, data type, and analysis method from the list below,</p>
            <p>to get hints on the calling sequence and parameters:</p>
        </div>
    """

    if is_notebook():
        try:
            # For scripts
            current_dir = os.path.dirname(os.path.abspath(__file__))
        except NameError:
            # For Jupyter notebooks
            current_dir = os.getcwd()
        img_path = os.path.join(current_dir, '..', 'assets', 'WaLSAtools_black.png')
        # display(HTML(f'<img src="{img_path}" style="margin-left: 40px; margin-top: 20px; width:300px; height: auto;">')) # not shwon in Jupyter notebook, only in MS Code
        import base64
        # Convert the image to Base64
        with open(img_path, "rb") as img_file:
            encoded_img = base64.b64encode(img_file.read()).decode('utf-8')
        # Embed the Base64 image in the HTML
        html_code_logo = f"""
        <div style="margin-left: 30px; margin-top: 20px;">
            <img src="data:image/png;base64,{encoded_img}" style="width: 300px; height: auto;">
        </div>
        """
        display(HTML(html_code_logo))
        display(HTML(credits_notebook))
    else:
        print(logo_terminal)
        print(credits_terminal)

from .parameter_definitions import display_parameters_text, single_series_parameters, cross_correlation_parameters

# Terminal-based interactive function
import textwrap
import shutil
def walsatools_terminal():
    """Main interactive function for terminal version of WaLSAtools."""
    print_logo_and_credits()

    # Get terminal width
    terminal_width = shutil.get_terminal_size().columns
    # Calculate 90% of the width (and ensure it's an integer)
    line_width = int(terminal_width * 0.90)

    # Step 1: Select Category
    while True:
        print("\n    Category:")
        print("    (a) Single time series analysis")
        print("    (b) Cross-correlation between two time series")
        category = input("    --- Select a category (a/b): ").strip().lower()

        if category not in ['a', 'b']:
            print("    Invalid selection. Please enter either 'a' or 'b'.\n")
            continue

        # Step 2: Data Type
        if category == 'a':
            while True:
                print("\n    Data Type:")
                print("    (1) 1D signal")
                print("    (2) 3D datacube")
                method = input("    --- Select a data type (1/2): ").strip()

                if method not in ['1', '2']:
                    print("    Invalid selection. Please enter either '1' or '2'.\n")
                    continue

                # Step 3: Analysis Method
                if method == '1':  # 1D Signal
                    while True:
                        print("\n    Analysis Method:")
                        print("    (1) FFT")
                        print("    (2) Wavelet")
                        print("    (3) Lomb-Scargle")
                        print("    (4) Welch")
                        print("    (5) EMD")
                        analysis_type = input("    --- Select an analysis method (1-5): ").strip()

                        method_map_1d = {
                            '1': 'fft',
                            '2': 'wavelet',
                            '3': 'lombscargle',
                            '4': 'welch',
                            '5': 'emd'
                        }

                        selected_method = method_map_1d.get(analysis_type, 'unknown')

                        if selected_method == 'unknown':
                            print("    Invalid selection. Please select a valid analysis method (1-5).\n")
                            continue

                        # Generate and display the calling sequence for 1D signal
                        print("\n    Calling sequence:\n")
                        return_values = single_series_parameters[selected_method]['return_values']
                        command = f">>> {return_values} = WaLSAtools(signal=INPUT_DATA, time=TIME_ARRAY, method='{selected_method}', **kwargs)"
                        wrapper = textwrap.TextWrapper(width=line_width - 4, initial_indent='    ', subsequent_indent='        ')
                        wrapped_command = wrapper.fill(command)
                        print(wrapped_command)

                        # Display parameter hints
                        display_parameters_text(selected_method, category='a')
                        return  # Exit the function after successful output

                elif method == '2':  # 3D Datacube
                    while True:
                        print("\n    Analysis Method:")
                        print("    (1) k-omega")
                        print("    (2) POD")
                        print("    (3) Dominant Freq / Mean Power Spectrum")
                        analysis_type = input("    --- Select an analysis method (1-3): ").strip()

                        if analysis_type == '3':  # Sub-method required
                            while True:
                                print("\n    Analysis method for Dominant Freq / Mean Power Spectrum:")
                                print("    (1) FFT")
                                print("    (2) Wavelet")
                                print("    (3) Lomb-Scargle")
                                print("    (4) Welch")
                                sub_method = input("    --- Select a method (1-4): ").strip()

                                sub_method_map = {'1': 'fft', '2': 'wavelet', '3': 'lombscargle', '4': 'welch'}
                                selected_method = sub_method_map.get(sub_method, 'unknown')

                                if selected_method == 'unknown':
                                    print("    Invalid selection. Please select a valid sub-method (1-4).\n")
                                    continue

                                # Generate and display the calling sequence for sub-method
                                print("\n    Calling sequence:\n")
                                return_values = 'dominant_frequency, mean_power, frequency, power_map'
                                command = f">>> {return_values} = WaLSAtools(data=INPUT_DATA, time=TIME_ARRAY, averagedpower=True, dominantfreq=True, method='{selected_method}', **kwargs)"
                                wrapper = textwrap.TextWrapper(width=line_width - 4, initial_indent='    ', subsequent_indent='        ')
                                wrapped_command = wrapper.fill(command)
                                print(wrapped_command)

                                # Display parameter hints
                                display_parameters_text(selected_method, category='a')
                                return  # Exit the function after successful output

                        method_map_3d = {
                            '1': 'k-omega',
                            '2': 'pod'
                        }

                        selected_method = method_map_3d.get(analysis_type, 'unknown')

                        if selected_method == 'unknown':
                            print("    Invalid selection. Please select a valid analysis method (1-3).\n")
                            continue

                        # Generate and display the calling sequence for k-omega/POD
                        print("\n    Calling sequence:\n")
                        return_values = single_series_parameters[selected_method]['return_values']
                        command = f">>> {return_values} = WaLSAtools(data1=INPUT_DATA, time=TIME_ARRAY, method='{selected_method}', **kwargs)"
                        wrapper = textwrap.TextWrapper(width=line_width - 4, initial_indent='    ', subsequent_indent='        ')
                        wrapped_command = wrapper.fill(command)
                        print(wrapped_command)

                        # Display parameter hints
                        display_parameters_text(selected_method, category='a')
                        return  # Exit the function after successful output

        elif category == 'b':  # Cross-correlation
            while True:
                print("\n    Data Type:")
                print("    (1) 1D signal")
                method = input("    --- Select a data type (1): ").strip()

                if method != '1':
                    print("    Invalid selection. Please enter '1'.\n")
                    continue

                while True:
                    print("\n    Analysis Method:")
                    print("    (1) FFT")
                    print("    (2) Wavelet")
                    print("    (3) Welch")
                    analysis_type = input("    --- Select an analysis method (1-2): ").strip()

                    cross_correlation_map = {
                        '1': 'fft',
                        '2': 'wavelet',
                        '3': 'welch'
                    }

                    selected_method = cross_correlation_map.get(analysis_type, 'unknown')

                    if selected_method == 'unknown':
                        print("    Invalid selection. Please select a valid analysis method (1-2).\n")
                        continue

                    # Generate and display the calling sequence for cross-correlation
                    print("\n    Calling sequence:\n")
                    return_values = cross_correlation_parameters[selected_method]['return_values']
                    command = f">>> {return_values} = WaLSAtools(data1=INPUT_DATA1, data2=INPUT_DATA2, time=TIME_ARRAY, method='{selected_method}', **kwargs)"
                    wrapper = textwrap.TextWrapper(width=line_width - 4, initial_indent='    ', subsequent_indent='        ')
                    wrapped_command = wrapper.fill(command)
                    print(wrapped_command)

                    # Display parameter hints
                    display_parameters_text(selected_method, category='b')
                    return  # Exit the function after successful output


# Jupyter-based interactive function
import ipywidgets as widgets # type: ignore
from IPython.display import display, clear_output, HTML # type: ignore
from .parameter_definitions import display_parameters_html, single_series_parameters, cross_correlation_parameters # type: ignore

# Global flag to prevent multiple observers
is_observer_attached = False

def walsatools_jupyter():
    """Main interactive function for Jupyter Notebook version of WaLSAtools."""
    global is_observer_attached, category, method, analysis_type, sub_method  # Declare global variables for reuse

    # Detach any existing observers and reset the flag
    try:
        detach_observers()
    except NameError:
        pass  # `detach_observers` hasn't been defined yet

    is_observer_attached = False  # Reset observer flag

    # Clear any previous output
    clear_output(wait=True)

    print_logo_and_credits()

    # Recreate widgets to reset state
    category = widgets.Dropdown(
        options=['Select Category', 'Single time series analysis', 'Cross-correlation between two time series'],
        value='Select Category',
        description='Category:'
    )
    method = widgets.Dropdown(
        options=['Select Data Type'],
        value='Select Data Type',
        description='Data Type:'
    )
    analysis_type = widgets.Dropdown(
        options=['Select Method'],
        value='Select Method',
        description='Method:'
    )
    sub_method = widgets.Dropdown(
        options=['Select Sub-method', 'FFT', 'Wavelet', 'Lomb-Scargle', 'Welch'],
        value='Select Sub-method',
        description='Sub-method:',
        layout=widgets.Layout(display='none')  # Initially hidden
    )

    # Persistent output widget
    output = widgets.Output()

    def clear_output_if_unselected(change=None):
        """Clear the output widget if any dropdown menu is unselected."""
        with output:
            if (
                category.value == 'Select Category'
                or method.value == 'Select Data Type'
                or analysis_type.value == 'Select Method'
                or (
                    analysis_type.value == 'Dominant Freq / Mean Power Spectrum'
                    and sub_method.layout.display == 'block'
                    and sub_method.value == 'Select Sub-method'
                )
            ):
                clear_output(wait=True)
                warn = '<div style="font-size: 1.1em; margin-left: 30px; margin-top:15px; margin-bottom: 15px;">Please select appropriate options from all dropdown menus.</div>'
                display(HTML(warn))

    def update_method_options(change=None):
        """Update available Method and Sub-method options."""
        clear_output_if_unselected()  # Ensure the output clears if any dropdown is unselected.
        sub_method.layout.display = 'none'
        sub_method.value = 'Select Sub-method'
        sub_method.options = ['Select Sub-method']

        if category.value == 'Single time series analysis':
            method.options = ['Select Data Type', '1D signal', '3D datacube']
            if method.value == '1D signal':
                analysis_type.options = ['Select Method', 'FFT', 'Wavelet', 'Lomb-Scargle', 'Welch', 'EMD']
            elif method.value == '3D datacube':
                analysis_type.options = ['Select Method', 'k-omega', 'POD', 'Dominant Freq / Mean Power Spectrum']
            else:
                analysis_type.options = ['Select Method']
        elif category.value == 'Cross-correlation between two time series':
            method.options = ['Select Data Type', '1D signal']
            if method.value == '1D signal':
                analysis_type.options = ['Select Method', 'FFT', 'Wavelet', 'Welch']
            else:
                analysis_type.options = ['Select Method']
        else:
            method.options = ['Select Data Type']
            analysis_type.options = ['Select Method']

    def update_sub_method_visibility(change=None):
        """Show or hide the Sub-method dropdown based on conditions."""
        clear_output_if_unselected()  # Ensure the output clears if any dropdown is unselected.
        if (
            category.value == 'Single time series analysis'
            and method.value == '3D datacube'
            and analysis_type.value == 'Dominant Freq / Mean Power Spectrum'
        ):
            sub_method.options = ['Select Sub-method', 'FFT', 'Wavelet', 'Lomb-Scargle', 'Welch']
            sub_method.layout.display = 'block'
        else:
            sub_method.options = ['Select Sub-method']
            sub_method.value = 'Select Sub-method'
            sub_method.layout.display = 'none'

    def update_calling_sequence(change=None):
        """Update the function calling sequence based on user's selection."""
        clear_output_if_unselected()  # Ensure the output clears if any dropdown is unselected.
        if (
            category.value == 'Select Category'
            or method.value == 'Select Data Type'
            or analysis_type.value == 'Select Method'
            or (
                analysis_type.value == 'Dominant Freq / Mean Power Spectrum'
                and sub_method.layout.display == 'block'
                and sub_method.value == 'Select Sub-method'
            )
        ):
            return  # Do nothing until all required fields are properly selected

        with output:
            clear_output(wait=True)

            # Handle Dominant Frequency / Mean Power Spectrum with Sub-method
            if (
                category.value == 'Single time series analysis'
                and method.value == '3D datacube'
                and analysis_type.value == 'Dominant Freq / Mean Power Spectrum'
            ):
                if sub_method.value == 'Select Sub-method':
                    print("Please select a Sub-method.")
                    return

                sub_method_map = {'FFT': 'fft', 'Wavelet': 'wavelet', 'Lomb-Scargle': 'lombscargle', 'Welch': 'welch'}
                selected_method = sub_method_map.get(sub_method.value, 'unknown')
                return_values = 'dominant_frequency, mean_power, frequency, power_map'
                command = f"{return_values} = WaLSAtools(signal=INPUT_DATA, time=TIME_ARRAY, averagedpower=True, dominantfreq=True, method='{selected_method}', **kwargs)"

            # Handle k-omega and POD
            elif (
                category.value == 'Single time series analysis'
                and method.value == '3D datacube'
                and analysis_type.value in ['k-omega', 'POD']
            ):
                method_map = {'k-omega': 'k-omega', 'POD': 'pod'}
                selected_method = method_map.get(analysis_type.value, 'unknown')
                parameter_definitions = single_series_parameters
                return_values = parameter_definitions.get(selected_method, {}).get('return_values', '')
                command = f"{return_values} = WaLSAtools(signal=INPUT_DATA, time=TIME_ARRAY, method='{selected_method}', **kwargs)"

            # Handle Cross-correlation
            elif category.value == 'Cross-correlation between two time series':
                cross_correlation_map = {'FFT': 'fft', 'Wavelet': 'wavelet', 'Welch': 'welch'}
                selected_method = cross_correlation_map.get(analysis_type.value, 'unknown')
                parameter_definitions = cross_correlation_parameters
                return_values = parameter_definitions.get(selected_method, {}).get('return_values', '')
                command = f"{return_values} = WaLSAtools(data1=INPUT_DATA1, data2=INPUT_DATA2, time=TIME_ARRAY, method='{selected_method}', **kwargs)"

            # Handle Single 1D signal analysis
            elif category.value == 'Single time series analysis' and method.value == '1D signal':
                method_map = {'FFT': 'fft', 'Wavelet': 'wavelet', 'Lomb-Scargle': 'lombscargle', 'Welch': 'welch', 'EMD': 'emd'}
                selected_method = method_map.get(analysis_type.value, 'unknown')
                parameter_definitions = single_series_parameters
                return_values = parameter_definitions.get(selected_method, {}).get('return_values', '')
                command = f"{return_values} = WaLSAtools(signal=INPUT_DATA, time=TIME_ARRAY, method='{selected_method}', **kwargs)"

            else:
                print("Invalid configuration.")
                return

            # Generate and display the command in HTML
            html_code = f"""
            <div style="font-size: 1.2em; margin-left: 30px; margin-top:15px; margin-bottom: 15px;">Calling sequence:</div>
            <div style="display: flex; margin-left: 30px; margin-bottom: 3ch;">
                <span style="color: #222; min-width: 4ch;">>>> </span>
                <pre style="
                    white-space: pre-wrap; 
                    word-wrap: break-word;  
                    color: #01016D; 
                    margin: 0;
                ">{command}</pre>
            </div>
            """
            display(HTML(html_code))
            display_parameters_html(selected_method, category=category.value)

    def attach_observers():
        """Attach observers to the dropdown widgets and ensure no duplicates."""
        global is_observer_attached
        if not is_observer_attached:
            detach_observers()
            category.observe(clear_output_if_unselected, names='value')
            method.observe(clear_output_if_unselected, names='value')
            analysis_type.observe(clear_output_if_unselected, names='value')
            category.observe(update_method_options, names='value')
            method.observe(update_method_options, names='value')
            analysis_type.observe(update_sub_method_visibility, names='value')
            sub_method.observe(update_calling_sequence, names='value')
            analysis_type.observe(update_calling_sequence, names='value')
            is_observer_attached = True

    def detach_observers():
        """Detach all observers to prevent multiple triggers."""
        try:
            category.unobserve(update_method_options, names='value')
            method.unobserve(update_method_options, names='value')
            analysis_type.unobserve(update_sub_method_visibility, names='value')
            sub_method.unobserve(update_calling_sequence, names='value')
            analysis_type.unobserve(update_calling_sequence, names='value')
        except ValueError:
            pass

    attach_observers()
    display(category, method, analysis_type, sub_method, output)


# Unified interactive function for both terminal and Jupyter environments
def interactive():
    if is_notebook():
        walsatools_jupyter()
    else:
        walsatools_terminal()