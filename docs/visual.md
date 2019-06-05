# Visualization in MEDYAN

## Syntax

### `list <list-object>`
Display the content of a list object

### `new <list-object>`
Append an item at the end of a list object

### `pause`
Pause the simulation thread at next checkpoint

### `quit`
Abort the program

### `resume`
Resume the paused simulation thread

### `set <attribute> <value>`
Set a variable to a certain value

### `show <attribute>`
Show the value of an attribute

## Possible user inputs
- `new visual.window`
- `list visual.window`
- `new visual.profile`
- `list visual.profile`
- `set visual.profile.0.target filament`
- `set visual.profile.3.target linker`
- `set visual.profile.2.color red`
- `del visual.profile.2`
- `set visual.profile.2.polygon.mode fill`
- `set visual.profile.0.path.mode tube`
- `set visual.profile.0.path.radius 5`
- `show visual.profile.0`
- `pause`
- `next`
- `resume`
- `quit`
