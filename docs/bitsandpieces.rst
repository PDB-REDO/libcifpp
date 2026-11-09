Bits & Pieces
=============

The *libcifpp* library offers some extra code that makes the life of developers a bit easier.

A progress bar
--------------

Applications based on *libcifpp* may have a longer run time. To give some feedback to the user running your application in a terminal you can use the :cpp:class:`cif::progress_bar`. This class will display an ASCII progress bar along with optional status messages, but only if output is to a real TTY (terminal).

A progress bar is also shown only if the duration is more than two seconds. To avoid having flashing progress bars for short actions.

The progress bar uses an internal progress counter that starts at zero and ends when the max value has been reached after which it will be removed from the screen. Updating this internal progress counter can be done by adding a number of steps calling :cpp:func:`cif::progress_bar::consumed` or by setting the exact value for the counter by calling :cpp:func:`cif::progress_bar::progress`.

Colouring output
----------------

It is also nice to emphasise some output in the terminal by using colours. For this you can create output manipulators using :cpp:func:`cif::coloured`. To write a string in white, and bold letters on a red background you can do:

.. code-block:: cpp

	using namespace cif::colour;
	std::cout << cif::coloured("Hello, world!", white, red, bold) << '\n';

