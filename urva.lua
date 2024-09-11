load("apptainer")

LmodMessage("Make sure that `@BASEPATH = \"/urva/\"` is set in your input file.")
LmodMessage("Usage is simply: `urva inp.inp`.")

local img_name      = 'urva_0.2.0_amd64_2023_09_29_11_18_06.sif'
local img_directory = '/work/group/catco/applications/urva'
local img_path      = pathJoin(img_directory, img_name)

function build_command(cmd)
  local cmd_beginning = 'apptainer exec '
  local cmd_ending    = img_path .. ' python3 /urva/main.py '
  local sh_ending     = '"$@"'
  local csh_ending    = '$*'
  local sh_cmd        = cmd_beginning .. cmd_ending .. sh_ending
  local csh_cmd       = cmd_beginning .. cmd_ending .. csh_ending
  set_shell_function(cmd, sh_cmd, csh_cmd)
end

build_command("urva")

