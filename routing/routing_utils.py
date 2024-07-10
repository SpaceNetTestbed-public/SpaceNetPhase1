def generate_ips_for_constellation(ip_subnet):
    """
    Generate a list of /28 networks starting at the given ip_subnet value.

    Args:
        ip_subnet (str): The IP subnet in the format "x.x.x.x".

    Returns:
        list: A list of available IP addresses in the format [(1, ip_address)]; '1' indicates the IP is available.
    """
    ip_octs = ip_subnet.split(".") #['192', '168', '100', '0']
    available_ips = [] #initialize list of available IPs
    for i in range (int(ip_octs[1]), 250):  # starting with the first IP in the second octet up to 249
        for j in range (int(ip_octs[2]), 250): # starting with the first IP in the third octet up to 249
            for k in range (0, 240, 16): # k starts at 0 and goes up to 240 in increments of 16
                ip = str(str(ip_octs[0])+".")+str(i)+"."+str(j)+"."+str(k) # create the IP address
                available_ips.append((1, ip)) # add the IP address to the list of available IPs

    return available_ips # return the list of available IPs

def generate_ips_for_physical_nodes(num_of_ips):
    """
    Generate a list of 'num_of_ips' 192.168.101.x IP addresses for physical nodes.

    Args:
        num_of_ips (int): The number of IP addresses to generate.

    Returns:
        list: A list of tuples containing the index and IP address.

    Example:
        >>> generate_ips_for_physical_nodes(5)
        [(1, '192.168.101.1'), (1, '192.168.101.2'), (1, '192.168.101.3'), (1, '192.168.101.4')]; '1' indicates the IP is available.
    """
    available_ips = [] #initialize list of available IPs
    for i in range(1, num_of_ips): 
        ip = str("192.168.101.")+str(i) # create the IP address
        available_ips.append((1, ip)) # add the IP address to the list of available IPs

    return available_ips # return the list of available IPs

def get_free_network_address(pool):
    """
    Get the first available network address from the given pool and removes that address from the pool.

    Args:
        pool (list): A list of tuples representing the network pool. Each tuple contains two elements:
                     the first element indicates the availability (1 for available),
                     and the second element is the network address.

    Returns:
        int: The first available network address from the pool. Returns -1 if no address is available.
    """
    free_ip = -1 #initialize free_ip to -1
    for i in pool: #iterate through the pool
        if i[0] == 1: #if the first element of the tuple is 1
            free_ip = i[1] #set free_ip to the second element of the tuple
            pool.remove(i) #remove the tuple from the pool
            break#; #break the loop

    return free_ip #return the free_ip

def get_network_address(str_ip_address):
    """
    Returns the network address for a given IP address assuming a /28 subnet mask.

    Parameters:
    str_ip_address (str): The IP address in string format.

    Returns:
    str: The network address in string format.
    """
    ip_oct1, ip_oct2, ip_oct3, ip_oct4 = str_ip_address.split(".") #split the IP address into octets
    net_add1= int(ip_oct1) & 255 #calculate the network address for the first octet
    net_add2= int(ip_oct2) & 255 #calculate the network address for the second octet
    net_add3= int(ip_oct3) & 255 #calculate the network address for the third octet
    net_add4= int(ip_oct4) & 240 #calculate the network address for the fourth octet

    return str(net_add1)+"."+str(net_add2)+"."+str(net_add3)+"."+str(net_add4) #return the network address

def assign_ips_for_constellation(links, addresses_pool):
    """
    Assigns IP addresses for each link in a constellation.

    Args:
        links (list): List of links in the constellation.
        addresses_pool (list): List of available IP addresses.

    Returns:
        list: List of dictionaries containing the interface number as the key and assigned IP address for the interface as the value. Each link has two dictionaries and use the first 2 IPs in the /28 network.
    """
    list_of_Intf_IPs = [] #initialize list of interface IPs
    #link = satx-ethy:satz-ethw
    for link in links: #iterate through the links
        linkIntf1, linkIntf2 = link.split(":") #split the link into the component interfaces
        network_address = get_free_network_address(addresses_pool) #get the first available network address from the pool
        if network_address != -1: #if a network address is available
           oct1, oct2, oct3, oct4 = network_address.split('.')#; #split the network address into octets
           list_of_Intf_IPs.append({"Interface": linkIntf1, "IP": oct1+"."+oct2+"."+oct3+"."+str(int(oct4)+1)+"/28"}) #add the interface and IP to the list of interface IPs
           list_of_Intf_IPs.append({"Interface": linkIntf2, "IP": oct1+"."+oct2+"."+oct3+"."+str(int(oct4)+2)+"/28"}) #add the interface and IP to the list of interface IPs
        else: #if no network address is available
           print("[Create Sat Network -- GSL] No Available IPs to assign") #print a message indicating no available IPs to assign
    return list_of_Intf_IPs #return the list of interface IPs

def get_link_intfs_ips(node1, node2, links, list_of_Intf_IPs):
    """
    UNUSED?
    Get the interface IPs of the links between two nodes.

    Args:
        node1 (str): Name of the first node.
        node2 (str): Name of the second node.
        links (list): List of links between nodes.
        list_of_Intf_IPs (list): List of interface IP dictionaries.

    Returns:
        list: List of interface IPs of the links between the two nodes.
    
    Note:
        Will error out if there is no link between the two given nodes
    """
    link_intfs_ips = [] #initialize list of interface IPs
    n1n2Link = "" #initialize n1n2Link to an empty string
    for link in links: #iterate through the links
        if node1 in link and node2 in link:  #if the nodes are in the link
            n1n2Link = link #set n1n2Link to the link
            break #break the loop

    print(n1n2Link) #print the link
    linkIntf1, linkIntf2 = n1n2Link.split(":") #split the link into the component interfaces
    # print linkIntf1, linkIntf2
    for intf_IP in list_of_Intf_IPs: #iterate through the list of interface IPs to find the IPs for the target interfaces
        if intf_IP["Interface"] == linkIntf1: #if the interface is the first interface of the link
            link_intfs_ips.append(intf_IP) #add the interface IP to the list of interface IPs
        if intf_IP["Interface"] == linkIntf2: #if the interface is the second interface of the link
            link_intfs_ips.append(intf_IP) #add the interface IP to the list of interface IPs

    return link_intfs_ips #return the list of interface IPs

def get_node_intf_ip(interface, list_of_Intf_IPs):
    """
    Get the IP address associated with the given interface from the list of interface IPs.

    Args:
        interface (str): The interface name.
        list_of_Intf_IPs (list): A list of dictionaries containing interface and IP information.

    Returns:
        str: The IP address associated with the given interface.

    Note:
        Will return nothing if the interface is not found in the list.
    """

    for intf_IP in list_of_Intf_IPs: #iterate through the list of interface IPs
        if intf_IP["Interface"] == interface: #if the interface is found
            return intf_IP["IP"] #return the IP address associated with the interface

    for intf_IP in list_of_Intf_IPs:
        if intf_IP["Interface"] == interface:
            return intf_IP["IP"]
        
def prepare_routing_config_commands(topology, data_path, initial_routes, links, list_of_Intf_IPs, satellites_by_index, num_of_threads):
    """
    This function prepares routing configuration commands by creating static routes in parallel, 
    logging the commands, and then writing each command to a separate file in a specific directory.

    Args:
        topology (object):          Sat_network class object. 
        data_path (str):            Path for simulation results
        initial_routes (dict):      Initial routes in the topology
        links (dict):               Corresponding links for the node pairs in the topology
        list_of_Intf_IPs (dict):    Assigned IPs for each interface
        satellites_by_index (dict): Dictionary of satellites by index
        num_of_threads (int):       Number of CPU threads for multi-threading

    Returns:
        None:                       Writes the commands to files but does not return anything
    """

    # Create static routes in parallel and store the commands in ipRouteCMD
    ipRouteCMD = topology.create_static_routes_batch_parallel(initial_routes, links, list_of_Intf_IPs, satellites_by_index, num_of_threads)
    
    # Open a log file to write the commands
    logg = open(data_path+"/stat_r.sh", "w")
    for rcmd in ipRouteCMD:
        for c in rcmd:
            # Write each command to the log file
            logg.write(c)

    # Close the log file
    logg.close()

    # Open the log file in read mode
    file1 = open(data_path+"/stat_r.sh", 'r')

    # Read all lines from the log file
    Lines = file1.readlines()

    # If the directory cmd_files does not exist, create it
    if os.path.isdir(data_path+"/cmd_files") == False:
        os.mkdir(data_path+"/cmd_files")

    # Remove all files in the cmd_files directory
    for f in os.listdir(data_path+"/cmd_files"):
        os.remove(os.path.join(data_path+"/cmd_files", f))

    # For each line in the log file
    for line in Lines:

        # Split the line into a list of commands
        command = line.strip().split(" ")

        # Open a new file in append mode for each command
        file = open(data_path+"/cmd_files/"+command[0]+"_routes.sh", 'a')
        string_to_write = ""

        # Concatenate all parts of the command except the first one
        for i in range(1,len(command)):
            string_to_write += command[i]+" "

        # Write the concatenated command to the file
        file.writelines(string_to_write+"\n")

        # Close the file
        file.close()
