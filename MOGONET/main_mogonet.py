from train_test import train_test

if __name__ == "__main__":    
    data_folder = './'
    view_list = [1,2,3]
    num_epoch_pretrain = 500
    num_epoch = 500
    lr_e_pretrain = 1e-3
    lr_e = 5e-4
    lr_c = 1e-3
    
    num_class = 2
    
    train_test(data_folder, view_list, num_class,
               lr_e_pretrain, lr_e, lr_c, 
               num_epoch_pretrain, num_epoch)             
