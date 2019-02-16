# Execute the example file
- Compile ROPTLIB as a library
    - ```make libropt.so```
- Copy libropt.so into PythonExamples folder
- Compile Cython
    - ```python setup2.py build_ext --inplace```
- Use in Python
    ```
    from PySimpleExample import py_simple_example
    result = py_simple_example()
     ```

# Create out your own code for using existing 
- These instructions assume that
    - you want to use one of the existing problems (cost function + gradient + option Hessian).
    - all  you need to do is to combine existing ROPTLIB objects in a particular combination of problem+solver+parameters and, perhaps, add some user defined code.
- Steps:
    - Make sure that all the classes and methods that you need are included in libropt.pxd. Use existing file as template.
    - Write your Cython wrapper in a .pyx file. You can ese PySimpleExample.pyx as a template. The wrapper will provide the connection between ROPTLIB C++ code and Python.
    - Compile your wrapper.
    - Use from Python.



