# cloudy_tools
Pythonic tools for working with cloudy. Simply download cloudy_outputs.py and place it in your cloudy workspace.

See example notebook for building cloudy models with custom stellar populations.

Basic useage:
- Place cloudy_outputs.py in your cloudy workspace, where output files are stored.
- Import cloudy_outputs into your workspace of choice.
- Use cloudy_outputs.get_cloudy_table('root'), where 'root' is the root name of all your cloudy outputs.
- This produces a table that can be manipulated with e.g astropy.
