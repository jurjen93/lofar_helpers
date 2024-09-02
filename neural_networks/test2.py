import torch
import torch.distributed as dist
import torch.multiprocessing as mp
import torch.nn as nn
import torch.optim as optim
import os
from torch.nn.parallel import DistributedDataParallel as DDP

def example(rank, local_rank, world_size):
    print(rank, local_rank, world_size)
    # Initialize process group
    torch.cuda.set_device(local_rank)
    
    dist.init_process_group("nccl", rank=rank, world_size=world_size)
    dist.barrier()
    
    # # Set the device for this process
    # Create model and move it to the appropriate GPU
    model = nn.Linear(10, 10).to(local_rank)
    
    # Construct DDP model
    ddp_model = DDP(model, device_ids=[local_rank])
    
    # Define loss function and optimizer
    loss_fn = nn.MSELoss()
    optimizer = optim.SGD(ddp_model.parameters(), lr=0.001)
    
    # Forward pass
    outputs = ddp_model(torch.randn(20, 10).to(local_rank))
    labels = torch.randn(20, 10).to(local_rank)
    
    # Backward pass
    loss = loss_fn(outputs, labels)
    loss.backward()
    
    # Synchronize gradients across all processes
    param_grad = next(model.parameters()).grad
    grads = [torch.zeros_like(param_grad) for _ in range(dist.get_world_size())]
    dist.all_gather(grads, param_grad)

    # Check if all gathered gradients are the same
    synced = all(torch.equal(grads[0], g) for g in grads)
    if synced:
        print(f"Process {rank}: Gradients are synchronized across processes.")
    else:
        print(f"Process {rank}: Gradients are NOT synchronized!")

    # Update parameters
    optimizer.step()
    
    # # Synchronize all processes
    dist.barrier()
    # print(f"Process {rank}: Model weights - {model.weight}")

    dist.destroy_process_group()

    

if __name__ == "__main__":
    import socket
    import subprocess
    # Get SLURM environment variables for distributed training

    os.environ['MASTER_ADDR'] = subprocess.getoutput("scontrol show hostnames $SLURM_NODELIST | head -n 1")

    # os.environ['MASTER_PORT'] = '29400'  # You can use any free port here

    # SLURM environment variables

    # Get the world size from the WORLD_SIZE variable or directly from SLURM:
    world_size = int(os.environ.get('WORLD_SIZE', os.environ.get('SLURM_NTASKS')))
    # Likewise for RANK and LOCAL_RANK:
    rank = int(os.environ.get('SLURM_PROCID'))
    gpus_per_node = torch.cuda.device_count()
    local_rank = rank - gpus_per_node * (rank //gpus_per_node)
    print(torch.__version__)
    print('local_rank', local_rank)
    print('rank', rank)
    print(world_size)
    print()
    example(rank, local_rank, world_size)

    
