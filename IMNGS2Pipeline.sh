#!/usr/bin/env bash



printf "\t╦╔╦╗╔╗╔╔═╗╔═╗  ___  ╔═╗┬┌─┐┌─┐┬  ┬┌┐┌┌─┐\n"
printf "\t║║║║║║║║ ╦╚═╗  ___| ╠═╝│├─┘├┤ │  ││││├┤ \n"
printf "\t╩╩ ╩╝╚╝╚═╝╚═╝ |___  ╩  ┴┴  └─┘┴─┘┴┘└┘└─┘\n"

# validateAction checks if the actions if one of the known actions
function validateAction {
    # Validate action
    case $1 in
        # Basic docker commands
        build)          printf "| Action: Build\n";;
        build-no-cache)          printf "| Action: Build without docker cache\n";;
        up)             printf "| Action: Up\n";;
        scale)          printf "| Action: Scale\n";;
        stop)           printf "| Action: Stop\n";;
        rm)             printf "| Action: Remove\n";;
        help)            printf "| Action: Help\n";  printOptions; exit 1;;

        # Everything else
        "")             printOptions; exit 1;;
        *)              printf "| Invalid action $action\n";exit 1;;
    esac
}


function printOptions {
    printf "Usage:\tcrc.sh -a build\n";
    printf "\tcrc.sh -a [options] -s [options]\n\n";
    printf "a[ction] options:\n"
    printf "\tbuild\t\tBuilds all images.\n";
    printf "\tbuild-no-cache\tSame like build but invalidates docker cache. Build is much slower but containers are fresh. Normally used when build does not work correctly.\n";
    printf "\tscale\t\tScale the workers.Used for upscaling and downscaling.\n";
    printf "\tup\t\t(Re)Creates, starts, and attaches to containers for a service. Unless they are already running, this command also starts any linked services.\n\t\t\tEquivalent to the docker-compose command (see official manual for further options). Note: it won't be able to build the containers, use the build command instead.\n";
    printf "\tstop\t\tStops running containers without removing them. They can be started again with start. Equivalent to the docker-compose command (see official manual for further options).\n";
    printf "\trm\t\tRemoves stopped service containers.\n";
    printf "\thelp\t\tPrints this message.\n\n";
}


action=''
scale='1'
target='dev'

# Parse the flags
while getopts 'a:s:t' flag; do
  case "${flag}" in
    a) action="${OPTARG}" ;;
    s) scale="${OPTARG}" ;;
    t) target="${OPTARG}" ;;
    *) error "Unexpected option ${flag}" ;;
  esac
done


validateAction $action

compose_file="docker-compose-IMNGS2Pipeline.yml"

echo ${compose_file}
case $action in
    build)
        docker-compose -f $compose_file build
        ;;

    build-no-cache)
        docker-compose -f $compose_file build --no-cache > build.logs
        ;;

    up)
        docker-compose -f $compose_file up -d > up.logs
        ;;

    scale)
        docker-compose -f $compose_file up -d --scale worker=$scale > up.logs
        ;;

    stop)
        docker-compose -f $compose_file stop
        ;;
    
    rm)
        docker-compose rm
        ;;
esac