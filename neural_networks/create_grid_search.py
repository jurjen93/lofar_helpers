with open("parameters.txt", "w") as f:
    i = 0
    for model_name in ["dinov2_vitl14_reg", "dinov2_vitb14_reg", "efficientnet_v2_l"]:
        for lr, batch_size in [(1e-4, 32), (1e-4, 64), (2e-4, 64), (5e-5, 32)]:
            for dropout_p in [0.1, 0.15, 0.2, 0.25]:
                for label_smoothing in [0.1, 0.2, 0.3]:
                    i += 1
                    cmd = f"{model_name} {lr} 1 {dropout_p} {batch_size} {label_smoothing} 0 0 16 16 {560 if 'dinov2' in model_name else 0} stack 0\n"
                    print(cmd)
                    f.writelines(cmd)
print(i)
